import numpy as np
import math
from utils.Reference_system import *
from scipy.integrate import solve_ivp
from integrator.Integrators import RK4
from integrator.Events import *
import time
from datetime import *
from utils.closed_form_solution import closed_form
import config
from control.heatload_control.utils_timeswitch import *


def asim(ip, m, time_0, OE, args, k_cf, heat_rate_control, time_switch_2 = 0,  reevaluation_mode = 1, time_switch_eval=False):
    # from physical_models.Gravity_models import gravity_invsquaredandJ2effect as gm
    if ip.gm == 0:
        from physical_models.Gravity_models import gravity_const as gm
    elif ip.gm == 1:
        from physical_models.Gravity_models import gravity_invsquared as gm
    elif ip.gm == 2:
        from physical_models.Gravity_models import gravity_invsquaredandJ2effect as gm

    # Aerodynamic Model
    from physical_models.Aerodynamic_models import aerodynamicscoefficient_fM as am

    # Thermal Model
    from physical_models.Thermal_models import heatrate_convective_maxwellian as hr_c

    version = args.MarsGram_version

    [r0, v0] = orbitalelemtorv(OE, m.planet)
    if config.count_numberofpassage != 1:
        t_prev = config.solution.orientation.time[-1]
    else:
        t_prev = m.initialcondition.time_rot
    v0_pp = r_intor_p(r0, v0, m.planet, time_0, t_prev)[1]
    date_initial = datetime(year=m.initialcondition.year,month=m.initialcondition.month,day=m.initialcondition.day,hour=m.initialcondition.hour,minute=m.initialcondition.min,second=m.initialcondition.second)

    T = m.planet.T  # fixed temperature
    RT = T * m.planet.R
    S = (np.linalg.norm(v0_pp) / (2 * RT) ** 0.5)
    CL_90, CD_90 = am(np.pi / 2, m.body, T, S, m.aerodynamics,montecarlo=0)  # Cl and Cd for all the body perpendicular
    CL_0, CD_0 = am(0, m.body, T, S, m.aerodynamics, montecarlo=0)  # Cl and Cd for all the body perpendicular
    CD_slope = (CD_90 - CD_0) / (math.pi * 0.5)

    def f(t0, in_cond, m):
        # Clock
        time_real = date_initial + timedelta(seconds=t0)
        timereal = clock(time_real.year, time_real.month, time_real.day, time_real.hour, time_real.minute, time_real.second)

        # Assign states
        pos_ii = in_cond[0:3]  # Inertial position
        pos_ii += 0.
        vel_ii = in_cond[3:6]  # Inertial velocity
        vel_ii += 0.
        mass = m.body.Mass  # Mass, kg
        pos_ii_mag = np.linalg.norm(pos_ii)  # Inertial position magnitude
        vel_ii_mag = np.linalg.norm(vel_ii)  # Inertial velocity magnitude
        lambdav_ii = in_cond[6]
        lambdagamma_ii = in_cond[7]
        lambdah_ii = in_cond[8]

        # Assign parameters
        omega_planet = m.planet.omega
        gamma = m.planet.gamma
        area_tot = m.body.Area_tot

        # TRANSFORM THE STATE
        # Inertial to planet relative transformation
        [pos_pp, vel_pp] = r_intor_p(pos_ii, vel_ii, m.planet, t0,t_prev) # Position vector planet / planet[m] and Velocity vector planet / planet[m / s]
        pos_pp_mag = np.linalg.norm(pos_pp)  # m, magnitude of position vector in PCPF

        vel_pp_mag = np.linalg.norm(vel_pp)

        # Orbital Elements
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)

        # Angular Momentum Calculations
        h_ii = np.cross(pos_ii, vel_ii)  # Inertial angular momentum vector[m ^ 2 / s]
        index = 0
        for item in h_ii:
            if item == -0:
                h_ii[index] = 0
                index = index + 1
        h_ii_mag = np.linalg.norm(h_ii)  # Magnitude of inertial angular momentum vector[m ^ 2 / s]
        h_pp = np.cross(pos_pp, vel_pp)

        index = 0
        for item in h_pp:
            if item == -0:
                h_pp[index] = 0
                index = index + 1
        h_pp_mag = np.linalg.norm(h_pp)
        h_pp_hat = h_pp / h_pp_mag

        # Inertial flight path angle
        arg = np.median([-1, 1, h_ii_mag / (pos_ii_mag * vel_ii_mag)])  # limit to[-1, 1]
        gamma_ii = math.acos(arg)
        if np.inner(pos_ii, vel_ii) < 0:
            gamma_ii = -gamma_ii

        # Relative flight path angle
        arg = np.median([-1, 1, h_pp_mag / (pos_pp_mag * vel_pp_mag)])  # limit to[-1, 1]
        gamma_pp = math.acos(arg)
        if np.inner(pos_pp, vel_pp) < 0:
            gamma_pp = -gamma_pp

        # Derived Quantity Calculations

        # Compute latitude and longitude
        LatLong = rtolatlong(pos_pp, m.planet)
        lat = LatLong[1]
        lon = LatLong[2]
        alt = LatLong[0]

        # Compute NED basis unit vectors
        uDuNuE = latlongtoNED(LatLong)  # nd
        uD = uDuNuE[0]
        uE = uDuNuE[2]
        uN = uDuNuE[1]

        # Get density, pressure, temperature and winds
        rho,T_p,wind = dm(h=alt, p=m.planet, OE=OE, lat=lat, lon=lon, timereal=timereal, t0=t0, tf_prev = t_prev, montecarlo=0,Wind=True, args=args, version = version)

        # Mach number
        sound_velocity = (gamma * m.planet.R * T_p) ** 0.5
        Mach = vel_pp_mag / sound_velocity
        S = ((gamma*0.5) ** 0.5) * Mach  # molecular speed ratio

        # CD = am(aoa=temp, body=m.body, T=T_p, S=S, args=args, montecarlo=0)[1]  # for the lambda_switch calculation, considering prev angle (!Approx!)

        # CD_90 = am(aoa=np.pi/2, body=m.body, T=T_p, S=S, args=args, montecarlo=0)[1]  # Cl and Cd for all the body perpendicular
        # CD_0 = am(aoa=0, body=m.body, T=T_p, S=S, args=args, montecarlo=0)[1]  # Cl and Cd for all the body perpendicular
        # CD_slope = (CD_90-CD_0)/(math.pi/2)


        if time_switch_eval == True:

            lambda_switch = np.divide(k_cf  * 2.0 * m.body.Mass * vel_ii_mag,area_tot * CD_slope * np.pi)
            if args.heat_load_sol == 0:
                if lambdav_ii < lambda_switch: #correct
                   aoa = 0.0001
                else:
                    aoa = m.aerodynamics.aoa
            elif args.heat_load_sol == 1:
                if lambdav_ii < lambda_switch:
                    aoa = m.aerodynamics.aoa
                else:
                    aoa = 0.0001
        else:
            if (args.heat_load_sol == 0 or args.heat_load_sol== 3):
                if t0 >= config.time_switch_1 and t0 <= time_switch_2: #correct
                    aoa = 0.0001
                else:
                    aoa = m.aerodynamics.aoa
            elif args.heat_load_sol == 1 or args.heat_load_sol== 2:
                if t0 >= config.time_switch_1 and t0 <= time_switch_2:
                    aoa = m.aerodynamics.aoa
                else:
                    aoa = 0.0001

        # Heat Rate
        heat_rate = hr_c(S, T_p, m, rho, vel_pp_mag,aoa)

        ## Add the control for the heat rate if flash == 3
        if heat_rate_control == True and heat_rate > args.max_heat_rate:
            from control.Control import control_solarpanels_heatrate
            state = [T_p, rho,S]
            index_ratio = [1]
            aoa = control_solarpanels_heatrate(ip,m,index_ratio,state)
            heat_rate = args.max_heat_rate

        # print('heat rate', heat_rate,'heat load', in_cond[-1])
        # print('heat load', in_cond[-1])
        # print('heat rate',heat_rate)

        # Convert wind to pp(PCPF) frame
        wE = wind[0] # positive to the east , m / s
        wN = wind[1] # positive to the north , m / s
        wU = wind[2] # positive up , m / s

        wind_pp = wN * uN + wE * uE - wU * uD  # wind velocity in pp frame , m / s
        vel_pp_rw = vel_pp + wind_pp  # relative wind vector , m / s
        vel_pp_rw_hat = vel_pp_rw / np.linalg.norm(vel_pp_rw)  # relative wind unit vector , nd
        # Dynamic pressure, CHANGE THE VELOCITY WITH THE WIND VELOCITY
        q = 0.5 * rho * np.linalg.norm(vel_pp_rw) ** 2  # base on wind - relative velocity

        ## Rotation Calculation
        rot_angle = np.linalg.norm(omega_planet) * t0  # rad
        L_PI = [[math.cos(rot_angle), math.sin(rot_angle), 0.0],
                [-math.sin(rot_angle), math.cos(rot_angle), 0.0],
                [0.0, 0.0, 1.0]]

        L_PI = [[+0. if x == -0 else x for x in row] for row in L_PI]

        g_ii = gm(pos_ii_mag, pos_ii, p=m.planet, mass=mass, vel_ii=vel_ii)
        gravity_ii = mass * g_ii


        bank_angle = 0.0
        lift_pp_hat = np.cross(h_pp_hat,vel_pp_rw_hat) #perpendicular vector to angular vector and velocity

        # Vehicle Aerodynamic Forces
        # CL and CD
        [CL,CD] = am(aoa=aoa, body=m.body, T=T_p, S=S, args=args, montecarlo=0)
        # Force calculations
        drag_pp_hat = -vel_pp_rw_hat # Planet relative drag force direction


        drag_pp = q * CD * area_tot * drag_pp_hat  # Planet relative drag force vector
        lift_pp = q * CL * area_tot * lift_pp_hat* math.cos(bank_angle)  # Planet relative lift force vector

        drag_ii = np.inner(np.transpose(L_PI), drag_pp)  # Inertial drag force vector
        lift_ii = np.inner(np.transpose(L_PI), lift_pp)  # Inertial lift force vector

        # Total Force
        # Total inertial external force vector on body [N]
        force_ii = drag_ii + lift_ii + gravity_ii

        # g_ii = np.divide(m.planet.g_ref * m.planet.Rp_e ** 2, np.square(pos_ii_mag))
        g_ii = np.linalg.norm(g_ii)
        # EOM

        # CD = CD_0+aoa*CD_slope
        lambdav_dot = -3 * k_cf * rho * vel_ii_mag ** 2 * aoa/ math.pi + lambdav_ii * (
                rho* area_tot * CD* vel_ii_mag) / mass - \
                      lambdagamma_ii * ((rho * area_tot * CL) / (2 * mass) + g_ii /  vel_ii_mag** 2 + 1 / (
                pos_ii_mag)) - lambdah_ii * gamma_ii
        lambdag_dot = lambdav_ii * g_ii - lambdah_ii * vel_ii_mag
        lambdah_dot = k_cf * rho * vel_ii_mag** 3 * aoa/ (math.pi * m.planet.H) - lambdav_ii * (
                (rho * area_tot *CD * vel_ii_mag ** 2) / (2 * mass * m.planet.H) + 2 * g_ii * gamma_ii/ (pos_ii_mag)) \
                      + lambdagamma_ii * (rho * area_tot * CL * vel_ii_mag / (2 * mass * m.planet.H) - 2 * g_ii / (
                (pos_ii_mag) * vel_ii_mag) + vel_ii_mag / (pos_ii_mag) ** 2)
        ydot = np.zeros(10)
        ydot[0:3] = vel_ii
        ydot[3:6] = force_ii / mass
        ydot[6] = lambdav_dot#-3*k_cf*rho*vel_ii_mag**2*aoa/math.pi + lambdav_ii* (rho * area_tot * CD *vel_ii_mag/ mass)-lambdagamma_ii* (rho * area_tot * CL/(2*mass) +g_ii/vel_ii_mag**2+1/pos_ii_mag)-lambdah_ii*gamma_ii
        ydot[7] = lambdag_dot#lambdav_ii*g_ii-lambdah_ii*vel_ii_mag
        ydot[8] = lambdah_dot#k_cf*rho*vel_ii_mag**3*aoa/(math.pi*m.planet.H)- \
                  #lambdav_ii* (rho * area_tot * CD *vel_ii_mag**2/ (2*mass*m.planet.H) + 2*g_ii*gamma_ii/pos_ii_mag)+\
                  #lambdagamma_ii*(rho*area_tot * CL*vel_ii_mag/(2*mass*m.planet.H)-2*g_ii/(pos_ii_mag*vel_ii_mag)+vel_ii_mag/(pos_ii_mag)**2)
        ydot[9] = heat_rate

        return ydot


    def out_drag_passage(t, y):
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - args.AE * 10 ** 3)
    out_drag_passage.terminal = True
    out_drag_passage.direction = 1

    def time_switch_fun(t, y):
        # time_real = date_initial + timedelta(seconds=t)
        # timereal = clock(time_real.year, time_real.month, time_real.day, time_real.hour, time_real.minute, time_real.second)
        # pos_ii = y[0:3]
        vel_ii = y[3:6]  # Inertial velocity
        vel_ii_mag = np.linalg.norm(vel_ii)  # Inertial velocity magnitude
        # [pos_pp,vel_pp] = r_intor_p(pos_ii, vel_ii, m.planet, t)
        # LatLong = rtolatlong(pos_pp, m.planet)
        # lat = LatLong[1]
        # lon = LatLong[2]
        # alt = LatLong[0]
        # rho,T_p,wind = dm(h=alt, p=m.planet, OE=OE, lat=lat, lon=lon, timereal=timereal, t0=t, montecarlo=0,Wind=True, args=args, version = version)

        # Mach number
        # sound_velocity = (m.planet.gamma * m.planet.R * T_p) ** 0.5
        # Mach = np.linalg.norm(vel_pp) / sound_velocity
        # S = ((m.planet.gamma*0.5) ** 0.5) * Mach  # molecular speed ratio
        # CD_90 = am(aoa=np.pi/2, body=m.body, T=T_p, S=S, args=args, montecarlo=0)[1]  # Cl and Cd for all the body perpendicular
        # CD_0 = am(aoa=0, body=m.body, T=T_p, S=S, args=args, montecarlo=0)[1]  # Cl and Cd for all the body perpendicular
        # CD_slope = (CD_90-CD_0)/(math.pi/2)
        lambda_switch = np.divide(k_cf * 2.0 * m.body.Mass * vel_ii_mag, m.body.Area_tot * CD_slope * np.pi)
        return (lambda_switch-y[6])
    time_switch_fun.terminal = False


    if time_switch_eval == True:
        # Density
        # from physical_models.Density_models import density_exp as dm

        if args.density_model == 'Exponential':
            from physical_models.Density_models import density_exp as dm
        elif args.density_model == 'MARSGram':
            from physical_models.Density_models import marsgram as dm

        # SOLVE EQUATIONS OF MOTIONS - 1 steps
        # USE CLOSED FORM SOLUTION TO DEFINE lambda_zero:
        T = m.planet.T  #fixed temperature
        t_cf, h_cf, gamma_cf, v_cf =closed_form(args,m,OE,T,aoa=m.aerodynamics.aoa,online=True)  # define closed-form solution
        #
        # RT = T * m.planet.R
        # S = (v_cf / (2 * RT) ** 0.5)
        # CL_90, CD_90 = am(np.pi / 2, m.body, T, S[0], m.aerodynamics,
        #                   montecarlo=0)  # Cl and Cd for all the body perpendicular
        # CL_0, CD_0 = am(0, m.body, T, S[0], m.aerodynamics, montecarlo=0)  # Cl and Cd for all the body perpendicular
        # CD_slope = (CD_90 - CD_0) / (math.pi * 0.5)
        #
        # coeff =  [CD_slope,CL_0, CD_0]
        # approx_sol = [t_cf, h_cf, gamma_cf, v_cf]

        # aoa_cf = aoa(m,k_cf, t_cf , h_cf , gamma_cf , v_cf,coeff)[0]
        # incond_lambda = func(k_cf, m, args, coeff, OE, heat_rate_control, approx_sol, aoa_cf, initial_guess=True)
        lambdav = v_cf[-1]#incond_lambda[0]#v_cf[-1]
        lambdag = 0.0#incond_lambda[1]#0.0
        lambdah = m.planet.mu / (m.planet.Rp_e + h_cf[-1]) ** 2#incond_lambda[2]#m.planet.mu / (m.planet.Rp_e + h_cf[-1]) ** 2

        lambda_v_fin = 10000
        lambda_h_fin = 10000
        lambda_gamma_fin = 10000

        lambda_v_fin_actual = 0
        lambda_gamma_fin_actual = 0
        lambda_h_fin_actual = 0
        count = 0
        while abs(lambda_v_fin_actual - lambda_v_fin) > 0.1 or abs(lambda_gamma_fin_actual - lambda_gamma_fin) > 0.1 or abs(lambda_h_fin_actual - lambda_h_fin) > 0.01:
            count += 1
            in_cond = [r0[0], r0[1], r0[2], v0[0], v0[1], v0[2], lambdav,lambdag,lambdah, 0]
            # Run Simulation
            eomfunction = lambda t, in_cond: f(t, in_cond, m)
            solution = solve_ivp(eomfunction, t_span=(time_0, time_0 + 1500), y0=in_cond, max_step=5,
                                 method='RK45',
                                 dense_output=True,
                                 events=[out_drag_passage])
            r_fin = (solution.y[0,-1]**2+solution.y[1,-1]**2+solution.y[2,-1]**2)**0.5
            v_fin = (solution.y[3,-1]**2+solution.y[4,-1]**2+solution.y[5,-1]**2)**0.5
            Q_fin = solution.y[-1,-1]
            lambda_v_fin = v_fin
            lambda_gamma_fin = 0.0
            lambda_h_fin = m.planet.mu / (r_fin)**2
            lambda_v_fin_actual =  solution.y[6,-1]
            lambda_gamma_fin_actual = solution.y[7,-1]
            lambda_h_fin_actual = solution.y[8,-1]
            in_cond = [solution.y[0,-1], solution.y[1,-1], solution.y[2,-1], solution.y[3,-1], solution.y[4,-1], solution.y[5,-1], lambda_v_fin,lambda_gamma_fin,lambda_h_fin, Q_fin]
            if (abs(lambda_v_fin_actual - lambda_v_fin) < 0.1 and abs(lambda_gamma_fin_actual - lambda_gamma_fin) < 0.1 and abs(lambda_h_fin_actual - lambda_h_fin) < 0.01) or count >4:
                break
            solution = solve_ivp(eomfunction, t_span=(solution.t_events[0],-10), y0=in_cond, max_step=5,
                                 method='RK45',
                                 dense_output=True,
                                 events=[out_drag_passage])
            # import matplotlib.pyplot as plt
            # plt.plot(solution.t, solution.y[6])
            # plt.show()

            lambdav = solution.y[6,-1]
            lambdag = solution.y[7,-1]
            lambdah = solution.y[8,-1]
            # print(abs(lambda_v_fin_actual - lambda_v_fin), abs(lambda_gamma_fin_actual - lambda_gamma_fin),
            #       abs(lambda_h_fin_actual - lambda_h_fin))
            # print('check',abs(lambda_v_fin_actual - lambda_v_fin) > 0.1 or abs(lambda_gamma_fin_actual - lambda_gamma_fin) > 0.1 or abs(lambda_h_fin_actual - lambda_h_fin) > 0.01)
            # print('check',abs(lambda_v_fin_actual - lambda_v_fin) > 0.1 and abs(lambda_gamma_fin_actual - lambda_gamma_fin) > 0.1 and abs(lambda_h_fin_actual - lambda_h_fin) > 0.01)

        # rerun the simulation with smaller step-size and the right lambda zero
        # Initial condition initialization
        in_cond = [r0[0], r0[1] , r0[2] , v0[0] , v0[1] , v0[2] , lambdav, lambdag, lambdah , 0]

        if args.density_model == 'Exponential':
            from physical_models.Density_models import density_exp as dm
        elif args.density_model == 'MARSGram':
            from physical_models.Density_models import marsgram as dm
        # from physical_models.Density_models import marsgram as dm
        eomfunction = lambda t , in_cond: f(t , in_cond , m)

        solution = solve_ivp(eomfunction , t_span=(time_0 , time_0 + 1500) , y0=in_cond , max_step=0.5,
                         method='RK45',
                         events=[out_drag_passage,time_switch_fun])
        temp = solution.t_events[1]


        # print('intersections', temp)
        ## time switch definition
        time_switch = [0, 0]
        if len(temp)==2:
            time_switch = temp
            # if abs(time_switch[1] - solution.t[-1]):
            #     time_switch[0] -= 0.1
        elif len(temp) == 1:
            time_switch[0] = temp
            time_switch[1] = solution.t[-1]


    else: # second time evaluation
        # Density
        from physical_models.Density_models import density_exp as dm1
        if args.density_model == 'Exponential':
            from physical_models.Density_models import density_exp as dm
        elif args.density_model == 'MARSGram':
            from physical_models.Density_models import marsgram as dm

        # SOLVE EQUATIONS OF MOTIONS - 1 steps


        temp_0 = 0
        tp = 1000
        # Initial condition initialization
        in_cond = [r0[0], r0[1] , r0[2] , v0[0] , v0[1] , v0[2], 0 ,0,0, config.heat_load_past]
        if reevaluation_mode == 1:  # bigger step for the first times revaluation is performed
            eomfunction = lambda t , in_cond: f(t , in_cond , m)
            solution = solve_ivp(eomfunction , t_span=(time_0 , time_0 + 1000) , y0=in_cond , max_step=1 ,
                                 method='RK23'  ,
                                 dense_output=False ,
                                 events=[out_drag_passage])
        elif reevaluation_mode == 2: # stricter conditions for the last times revaluation is performed
            eomfunction = lambda t , in_cond: f(t , in_cond , m)
            solution = solve_ivp(eomfunction , t_span=(time_0 , time_0 + 1000) , y0=in_cond , max_step=1 ,
                                 method='RK23',dense_output=False ,
                                 events=[out_drag_passage])
        return solution.y


    return solution.y, time_switch
