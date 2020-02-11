import numpy as np
import math
from utils.Reference_system import *
from scipy.integrate import solve_ivp
from integrator.Integrators import RK4
import config
from utils.save_results import *
from integrator.Events import *
import time
from datetime import *


def asim(ip, m, time_0, terminal_state, OE, numberofpassage, simulation):
    # Definition Models:
    # Gravity
    if ip.gm == 0:
        from physical_models.Gravity_models import gravity_const as gm
    elif ip.gm == 1:
        from physical_models.Gravity_models import gravity_invsquared as gm
    elif ip.gm == 2:
        from physical_models.Gravity_models import gravity_invsquaredandJ2effect as gm

    # Density
    if ip.dm == 0:
        from physical_models.Density_models import density_constant as dm
    elif ip.dm == 1:
        from physical_models.Density_models import density_exp as dm
    elif ip.dm == 2:
        from physical_models.Density_models import density_no as dm
    elif ip.dm == 3:
        from physical_models.Density_models import marsgram as dm

    # Wind
    wind_m = False
    if ip.wm == 1:
        wind_m = True

    # Aerodynamic Models
    if ip.am == 0:
        from physical_models.Aerodynamic_models import aerodynamicscoefficient_constant as am
    elif ip.am == 1:
        from physical_models.Aerodynamic_models import aerodynamicscoefficient_fM as am
    elif ip.am == 2:
        from physical_models.Aerodynamic_models import aerodynamicscoefficient_noballisticflight as am

    # Thermal Model
    if ip.tm.c == 1:
        from physical_models.Thermal_models import heatrate_convective as hr_c
    if ip.tm.r == 1:
        from physical_models.Thermal_models import heatrate_radiative as hr_r
    if ip.tm.t == 1:
        from physical_models.Thermal_models import heatrate_convective_maxwellian as hr_c

    # Control Model
    if ip.cm == 1:
        from physical_models.Control import control_solarpanels_openloop as cntr
    elif ip.cm == 0:
        from physical_models.Control import nocontrol as cntr

    # Monte Carlo Analysis
    MonteCarlo = False
    if ip.mc == 1:
        MonteCarlo = True

    if (OE.a > (3400 + 50 + 500) * 10 ** 3) and (simulation['Only DragPassage'] == False):
        index_steps_EOM = 3
    else:
        index_steps_EOM = 1
        r_p = OE.a * (1 - OE.e)
        print('Integration step number =', index_steps_EOM)

    # Rotation Matrix ijk tp pqw
    i = OE.i
    OMEGA = OE.OMEGA
    omega = OE.omega
    # matrix T_ijk: rotational matrix between RTN and PCI (satellite normal coordinate system)
    T_ijk = [[math.cos(OMEGA) * math.cos(omega) - math.sin(OMEGA) * math.sin(omega) * math.cos(i),
              math.sin(OMEGA) * math.cos(omega) + math.cos(OMEGA) * math.sin(omega) * math.cos(i),
              math.sin(omega) * math.sin(i)],
             [-math.cos(OMEGA) * math.sin(omega) - math.sin(OMEGA) * math.cos(omega) * math.cos(i),
              -math.sin(OMEGA) * math.sin(omega) + math.cos(OMEGA) * math.cos(omega) * math.cos(i),
              math.cos(omega) * math.sin(i)],
             [math.sin(OMEGA) * math.sin(i), -math.cos(OMEGA) * math.sin(i), math.cos(i)]]

    T_ijk = [[+0. if x == -0 else x for x in row] for row in T_ijk]
    [r0, v0] = orbitalelemtorv(OE, m.planet)

    config.count_numberofpassage = config.count_numberofpassage + 1
    def f(t0, in_cond, m, index_phase_aerobraking):
        ## Counters
        # Counter for all along the simulation of all passages
        config.count_aerobraking = config.count_aerobraking + 1
        passage_number = config.count_aerobraking
        # Counter for one entire passage
        config.count_dori = config.count_dori + 1
        # Counter for one phase
        config.count_phase = config.count_phase + 1


        # Clock

        time_real = datetime(year=m.initialcondition.year,month=m.initialcondition.month,day=m.initialcondition.day,hour=m.initialcondition.hour,minute=m.initialcondition.min,second=m.initialcondition.second) + timedelta(seconds=t0)
        timereal = clock(time_real.year, time_real.month, time_real.day, time_real.hour, time_real.minute, time_real.second)

        # Assign states
        pos_ii = in_cond[0:3]  # Inertial position
        pos_ii += 0.
        vel_ii = in_cond[3:6]  # Inertial velocity
        vel_ii += 0.
        r_i = cartesian(pos_ii[0], pos_ii[1], pos_ii[2])
        v_i = cartesian(vel_ii[0], vel_ii[1], vel_ii[2])
        mass = m.initialcondition.m  # Mass, kg
        pos_ii_mag = np.linalg.norm(pos_ii)  # Inertial position magnitude
        vel_ii_mag = np.linalg.norm(vel_ii)  # Inertial velocity magnitude

        # Assign parameters
        omega_planet = m.planet.omega
        gamma = m.planet.gamma
        mu_fluid = m.planet.mu_fluid
        area_tot = m.body.Area_SC + m.body.Area_SA

        # TRANSFORM THE STATE
        # Inertial to planet relative transformation
        [r_p, v_p] = r_intor_p(r_i, v_i, m.planet, t0)
        pos_pp = [r_p.x, r_p.y, r_p.z]  # Position vector planet / planet[m]
        pos_pp_mag = np.linalg.norm(pos_pp)  # m, magnitude of position vector in PCPF
        pos_pp_hat = pos_pp / pos_pp_mag  # nd, unit position vector in PCPF
        pos_ii_hat = pos_ii / pos_ii_mag  # Inertial position vector direction

        vel_pp = [v_p.x, v_p.y, v_p.z]  # Velocity vector planet / planet[m / s]
        vel_pp_mag = np.linalg.norm(vel_pp)
        vel_pp_hat = vel_pp / vel_pp_mag

        # Orbital Elements
        [OE] = rvtoorbitalelement(r_i, v_i, m.planet)

        if (OE.vi > 0 and OE.vi < math.pi/2) and config.ascending_phase == False:
            config.ascending_phase = True
            config.atmospheric_data = {}

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

        ## Derived Quantity Calculations

        # Compute latitude and longitude
        [LatLong] = rtolatlong(r_p, m.planet)
        lat = LatLong.LAT
        lon = LatLong.LON
        alt = LatLong.h

        if alt <= 160 * 10**3:
            config.drag_state = True
        else:
            config.drag_state = False

        # Compute NED basis unit vectors
        uDuNuE = latlongtoNED(LatLong)  # nd
        uD = uDuNuE.uD
        uE = uDuNuE.uE
        uN = uDuNuE.uN

        # compute azimuth
        vN = np.inner(vel_pp, uN)  # m / s
        vE = np.inner(vel_pp, uE)  # m / s
        azi_pp = math.atan2(vE, vN)  # rad

        # Get density, pressure, temperature and winds
        rho,T_p,wind = dm(OE, alt, lat, lon, timereal, t0, m.planet, montecarlo=MonteCarlo, Wind=wind_m)


        # Define output.txt containing density data
        #if (drag_state == True) and ():
        p = 0
        length_car = m.body.length_SA + m.body.length_SC
        Re = vel_pp_mag * rho * length_car / mu_fluid  # Reynolds number
        # Mach number
        a = (gamma * m.planet.R * T_p) ** 0.5
        M = vel_pp_mag / a
        S = ((gamma / 2) ** 0.5) * M  # molecular speed ratio

        if config.drag_state == True:
            ## Check type of fluid
            Kn = 1.26 * (gamma) ** 0.5 * M / Re
            if index_phase_aerobraking == 2 | 0:
                if (alt < 80000) & (config.index_warning_alt == 0):
                    print('WARNING: Altitude < 80 km!')
                    config.index_warning_alt = 1
                elif alt > 100000:
                    config.index_warning_alt = 0
                if (Kn < 0.1) & (config.index_warning_flow == 0):
                    print('Warning: transitional flow passage!')
                    config.index_warning_flow = 1
                elif Kn >= 0.1:
                    config.index_warning_flow = 0

        aoa = m.aerodynamics.aoa
        # Heat rate and Control
        if config.drag_state == True:
            # Theta definition
            #T_p = m.planet.T
            T_w = T_p
            #T_w = config.T_w
            aoa = cntr(m, rho, T_p, T_w, S, terminal_state)

            # Heat Rate
            T_r, St = hr_c(S, T_p, m.aerodynamics, m.planet)
            cp = m.planet.gamma / (m.planet.gamma - 1) * m.planet.R
            # heat_rate = ((m.aerodynamics.thermal_accomodation_factor*(m.planet.gamma + 1) /  m.planet.gamma) * (St*rho * vel_pp_mag * cp * (T_r - T_p)) )* 10**-4 # W/cm^2

            heat_rate = (m.aerodynamics.thermal_accomodation_factor * rho * m.planet.R * T_p) * (
                        (m.planet.R * T_p / (2 * math.pi)) ** 0.5) * (
                                    (S ** 2 + (gamma) / (gamma - 1) - (gamma + 1) / (2 * (gamma - 1)) * (T_w / T_p)) * (
                                        math.exp(-(S * math.sin(aoa)) ** 2) + (math.pi ** 0.5) * (S * math.sin(aoa)) *
                                        (1 + math.erf(S * math.sin(aoa)))) - 1 / 2 * math.exp(
                                -(S * math.sin(aoa)) ** 2)) * 10 ** -4  # W/cm^2
            #emissivity = 0.82
            #stef_bol_const = 5.67*10**-4
            #config.T_w = (heat_rate / (emissivity*stef_bol_const))**-4 - 4
            #print(config.T_w)
        else:
            T_r = 0
            heat_rate = 0

        # Convert wind to pp(PCPF) frame
        wE = wind[0] # positive to the east , m / s
        wN = wind[1] # positive to the north , m / s
        wU = wind[2] # positive up , m / s

        wind_pp = wN * uN + wE * uE - wU * uD  # wind velocity in pp frame , m / s
        vel_pp_rw = vel_pp + wind_pp  # relative wind vector , m / s
        vel_pp_rw_mag = np.linalg.norm(vel_pp_rw)  # relative wind magnitude , m / s
        vel_pp_rw_hat = vel_pp_rw / vel_pp_rw_mag  # relative wind unit vector , nd

        # Dynamic pressure, CHANGE THE VELOCITY WITH THE WIND VELOCITY
        q = 1 / 2 * rho * vel_pp_rw_mag ** 2  # base on wind - relative velocity

        ## Rotation Calculation
        rot_angle = np.linalg.norm(omega_planet) * t0  # rad
        L_PI = [[math.cos(rot_angle), math.sin(rot_angle), 0],
                [-math.sin(rot_angle), math.cos(rot_angle), 0],
                [0, 0, 1]]

        gravity_ii = mass * gm(pos_ii_mag, pos_ii, vel_ii, m.planet)


        # Vectors and Rotation
        # Tensors of interest
        n1 = np.array([[0, h_pp_hat[2], - h_pp_hat[1]],
                       [-h_pp_hat[2], 0, h_pp_hat[0]],
                       [h_pp_hat[1], -h_pp_hat[0], 0]])

        R1 = np.eye(3) + math.sin(math.pi/2) * n1 + (1 - math.cos(math.pi/2)) * np.inner(n1, n1)

        n2 = np.array([[0, vel_pp_rw_hat[2], - vel_pp_rw_hat[1]],
                       [-vel_pp_rw_hat[2], 0, vel_pp_rw_hat[0]],
                       [vel_pp_rw_hat[1], -vel_pp_rw_hat[0], 0]])  # change with vel_pp_rw_hat
        bank_angle =0
        R2 = np.eye(3) + math.sin(bank_angle) * n2 + (1 - math.cos(bank_angle)) * np.inner(n2, n2)
        # Vehicle Aerodynamic Forces
        # CL and CD
        [CL,CD] = am(T_p, S, m.aerodynamics, m.body, vel_pp_mag, rho, gamma_pp, aoa, montecarlo=MonteCarlo)

        # Force calculations
        drag_pp_hat = -vel_pp_rw_hat#vel_pp_rw_hat  # Planet relative drag force direction

        drag_pp = q * CD * area_tot * drag_pp_hat  # Planet relative drag force vector

        def perpendicular_vector(v):
            if v[1] == 0 and v[2] == 0:
                if v[0] == 0:
                    raise ValueError('zero vector')
                else:
                    return np.cross(v, [0, 1, 0])
            return np.cross(v, [1, 0, 0])

        lift_pp_hat = perpendicular_vector(vel_pp_rw_hat)  # Planet relative lift force direction
        #lift_pp_hat = np.inner(R2, np.inner(R1,vel_pp_rw_hat)) if reconsidered bank angle check before
        lift_pp = q * CL * area_tot * lift_pp_hat  # Planet relative lift force vector
        drag_ii = np.inner(np.transpose(L_PI), drag_pp)  # Inertial drag force vector

        lift_ii = np.inner(np.transpose(L_PI), lift_pp)  # Inertial lift force vector

        # Total Force
        # Total inertial external force vector on body [N]
        force_ii = drag_ii + lift_ii + gravity_ii

        # EOM
        ydot = np.zeros(6)

        ydot[0:3] = vel_ii
        ydot[3:6] = force_ii / mass
        energy = (vel_ii_mag ** 2) / 2 - m.planet.mu / (pos_ii_mag)

        # ## SAVE RESULTS
        sol = np.hstack([[t0], [timereal.year], [timereal.month], [timereal.day], [timereal.hour],[timereal.minute],
                         [timereal.second], [numberofpassage], pos_ii, vel_ii, [pos_ii_mag], [vel_ii_mag], pos_pp,
                         [pos_pp_mag], vel_pp, [vel_pp_mag], [OE.a], [OE.e], [OE.i], [OE.OMEGA], [OE.omega], [OE.vi],
                         [lat], [lon], [alt], [gamma_ii], [gamma_pp], h_ii, h_pp, [h_ii_mag], [h_pp_mag], uD, uE, uN,
                         [vN], [vE], [azi_pp], [rho], [T_p], [p], wind, [CL], [CD], [aoa], [mass], [heat_rate], [T_r], [q],
                         gravity_ii, drag_pp, drag_ii, lift_pp, lift_ii, force_ii, [energy], config.index_MonteCarlo, int(config.drag_state)]).reshape((1, 89)).tolist()
        config.solution_intermediate.extend(sol)
        return ydot

    ## EVENTS

    def eventfirststep(t, y):
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - 159 * 10 ** 3)

    eventfirststep.terminal = True

    def eventsecondstep(t, y):
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - 160 * 10 ** 3)

    eventsecondstep.terminal = True
    eventsecondstep.direction = 1

    def apoasispoint(t, y):
        y_pwj = np.inner(T_ijk, y[3:6])
        return y_pwj[0]

    apoasispoint.terminal = True
    apoasispoint.direction = 1

    def periapsispoint(t, y):
        y_pwj = np.inner(T_ijk, y[3:6])
        return y_pwj[0]

    periapsispoint.terminal = False
    periapsispoint.direction = -1

    def impact(t,
               y):  # nota: Event impact has to be the last because of the impact condition in the aerobraking script.
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - (
                    m.planet.Rp_e + 70000))  # consider impact when the altitude is less than 70 km

    impact.terminal = True

    def apoapsisgreaterperiapsis(t, y):
        r = y[0:3]
        v = y[3:6]

        Energy = (np.inner(v, v)) / 2 - m.planet.mu / (np.inner(r, r)) ** 0.5
        a = - m.planet.mu / (2 * Energy)
        h = np.cross(r, v)
        h += 0.
        e = (1 + (2 * Energy * (np.inner(h, h)) / (m.planet.mu) ** 2)) ** 0.5

        r_a = a * (1 + e)
        r_p = a * (1 - e)
        if r_a < (r_p):
            print('Periapsis greater than apoapsis!')
        return (r_a - r_p)

    apoapsisgreaterperiapsis.terminal = True

    #inizialize config parameter for new passage
    config.count_dori = 0
    config.atmospheric_data = {}
    config.ascending_phase = False
    config.state_inner_boundary_atmosphere = []


    ## SOLVE EQUATIONS OF MOTIONS - 1 steps
    if index_steps_EOM == 1: #(Spacecraft really close to the planet surface (semimajor axis small))
        # Initialization
        if time_0 == 0:
            clean_results()

        save_pre_index = len(config.solution.orientation.time)  # index definition pre drag passage
        index_phase_aerobraking = 0

        # Initial condition initialization
        in_cond = [r0.x, r0.y, r0.z, v0.x, v0.y, v0.z]
        y_0 = np.array(in_cond)
        if simulation['Only DragPassage'] == True:
            terminal_event = eventsecondstep
        else:
            terminal_event = apoasispoint

        # CONTROL IF
        if ip.cm == 0: # No control
            initial_time, final_time = time_0, 1e20
            t = [initial_time, final_time]

            # Run Simulation
            eomfunction = lambda t, in_cond: f(t, in_cond, m, index_phase_aerobraking)
            solution = solve_ivp(eomfunction, t, in_cond, method='RK45', dense_output=True, max_step=2, rtol=1e-8,
                                 atol=1e-10, events=[terminal_event, periapsispoint,
                                                     apoapsisgreaterperiapsis, impact])
            # Save results
            save_results(solution.t)

            # Check if hold requirements continue campaign
            continue_campaign = event(solution)

        elif ip.cm == 1: # Yes Control
            continue_simulation = True  # initialization breaker

            # Time initialization
            t_0, h, time_solution = time_0, 1, []

            # Run Simulation
            while continue_simulation:
                y, t, continue_simulation, continue_campaign = RK4(f, h, t_0, y_0, m, T_ijk, index_phase_aerobraking, simulation)

                # New initial condition
                y_0 = y
                t_0 = t

                # Save results
                time_solution.append(t)
            save_results(time_solution)

        save_post_index = len(config.solution.orientation.time)

        # Re-Set count index to 0
        config.count_dori = 0
        config.count_phase = 0

    ## SOLVE EQUATIONS OF MOTIONS - 3 steps
    # First Step
    else:
        # Initialization
        max_step_bigsize = 100
        index_phase_aerobraking = 1

        # Initial condition initialization
        in_cond = [r0.x, r0.y, r0.z, v0.x, v0.y, v0.z]
        # Time initialization
        initial_time, final_time = time_0, 1e20,
        t = [initial_time, final_time]

        # Run Simulation
        eomfunction = lambda t, in_cond: f(t, in_cond, m, index_phase_aerobraking)
        solution = solve_ivp(eomfunction, t, in_cond, method='RK45', dense_output=True, max_step=max_step_bigsize,
                             rtol=1e-5, atol=1e-8, events=[eventfirststep, periapsispoint, impact])

        # Save Results
        save_results(solution.t)
        save_pre_index = len(config.solution.orientation.time)

        # Re-Set count index to 0
        config.count_phase = 0
        print('End First step')

        # SECOND STEP
        # Initialization
        index_phase_aerobraking = 2

        # Initial condition initialization
        in_cond = [solution.y[0][-1], solution.y[1][-1], solution.y[2][-1], solution.y[3][-1], solution.y[4][-1],
                   solution.y[5][-1]]

        if ip.cm == 0:
            # Time initialization
            initial_time, final_time = solution.t[-1], 1e20
            t = [initial_time, final_time]


            # Run Simulation
            eomfunction = lambda t, in_cond: f(t, in_cond, m, index_phase_aerobraking)
            solution = solve_ivp(eomfunction, t, in_cond, method='RK45', dense_output=True, max_step=2, rtol=1e-8,
                                 atol=1e-10, events=[eventsecondstep, periapsispoint, impact])

            # Save Results
            save_results(solution.t)



        elif ip.cm == 1:
            continue_simulation = True  # initialization breaker

            y_0 = np.array(in_cond)  # Initial condition initialization

            # Time initialization
            t_0, h, time_solution = solution.t[-1], 1, []

            # Run Simulation
            while continue_simulation:
                y, t, continue_simulation, continue_campaign = RK4(f, h, t_0, y_0, m, T_ijk, index_phase_aerobraking, simulation)

                # New initial condition
                y_0 = y
                t_0 = t

                # Save results
                time_solution.append(t)
            save_results(time_solution)
        save_post_index = len(config.solution.orientation.time)

        # Re-Set count index to 0
        config.count_phase = 0
        print('End Second step')

        # THIRD STEP
        # Initialization
        index_phase_aerobraking = 3

        in_cond = [config.solution.orientation.pos_ii[0][-1], config.solution.orientation.pos_ii[1][-1],
                   config.solution.orientation.pos_ii[2][-1],
                   config.solution.orientation.vel_ii[0][-1], config.solution.orientation.vel_ii[1][-1],
                   config.solution.orientation.vel_ii[2][-1]]
        t = [config.solution.orientation.time[-1], final_time]

        # RUN SIMULATION
        eomfunction = lambda t, in_cond: f(t, in_cond, m, index_phase_aerobraking)
        solution = solve_ivp(eomfunction, t, in_cond, method='RK45', dense_output=True, max_step=max_step_bigsize,
                             rtol=1e-6, atol=1e-9, events=[apoasispoint, periapsispoint, impact])

        # Save results
        save_results(solution.t)

        # Re-Set count index to 0
        config.count_phase = 0
        print('End Third step')

        # Define breaker campaign
        continue_campaign = event(solution)


    ## Heat Load calculation (Because now time is independent by the simulation)
    heat_load = [None] * (len(config.solution.orientation.time) - config.save_index_heat)
    count = 0
    for counter in range(len(heat_load)):
        if config.solution.performance.heat_rate[counter + config.save_index_heat] > (m.aerodynamics.thermal_limit + 0.001):
            count += 1
        if (counter == 0):
            heat_load[counter] = 0
        elif config.solution.performance.heat_rate[counter + config.save_index_heat] == 0:
            heat_load[counter] = 0
        else:
            heat_load[counter] = heat_load[counter - 1] + (
                        config.solution.orientation.time[counter + config.save_index_heat] -
                        config.solution.orientation.time[counter - 1 + config.save_index_heat]) * \
                                 config.solution.performance.heat_rate[counter + config.save_index_heat]
    config.count_overcome_hr += count


    config.save_index_heat = len(config.solution.orientation.time)
    config.solution.performance.heat_load.extend(heat_load)
    ## Save Results # def new function
    print("Thermal limit overcomed {} times!".format(count))
    print("Actual periapsis altitude {0:.2f}km - Vacuum periapsis altitude = {1:.2f}km".format(min(config.solution.orientation.alt[save_pre_index:save_post_index]) / 10 ** 3, (config.solution.orientation.oe[0][-1] * (1 - config.solution.orientation.oe[1][-1]) - m.planet.Rp_e) / 10 ** 3))
    print("Ra new = {0:.2f}km".format((config.solution.orientation.oe[0][-1] * (1 + config.solution.orientation.oe[1][-1])) / 10 ** 3))
    config.altitudeperiapsis.append(min(config.solution.orientation.alt[save_pre_index:save_post_index]) / 10 ** 3)
    config.max_heatrate.append(max(config.solution.performance.heat_rate[save_pre_index:save_post_index]))

    return continue_campaign
