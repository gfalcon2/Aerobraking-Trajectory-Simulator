#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import numpy as np
import math
from utils.Reference_system import *
from scipy.integrate import solve_ivp
from integrator.Integrators import RK4
import config
from utils.save_results import *
from integrator.Events import *
from datetime import *
import numpy as np
cimport numpy as np
import cython
cimport cython
import traceback

from utils.Odyssey_maneuver_plan import Odyssey_firing_plan

def asim(ip, m, initial_state, numberofpassage, args):
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
    cdef int wind_m = False
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
    if ip.tm == 1:
        from physical_models.Thermal_models import heatrate_convective_radiative as hr_c
    if ip.tm == 2:
        from physical_models.Thermal_models import heatrate_convective_maxwellian as hr_c

    # Control Model
    if ip.cm == 3:
        from control.Control import control_solarpanels_openloop as cntr
    elif ip.cm == 2:
        from control.Control import control_solarpanels_heatload as cntr
    elif ip.cm == 1:
        from control.Control import control_solarpanels_heatrate as cntr
    elif ip.cm == 0:
        from control.Control import nocontrol as cntr

    # Thrust Maneuvers Model
    if ip.tc == 0:
        from control.Propulsive_maneuvers import no_manuever as thrust_m
    elif ip.tc == 1:
        from control.Propulsive_maneuvers import abms as thrust_m
    elif ip.tc == 2:
        from control.Propulsive_maneuvers import deceleration_drag_passage as thrust_m

    # Monte Carlo Analysis
    cdef int MonteCarlo = False
    if ip.mc == 1:
        MonteCarlo = True


    cdef int index_steps_EOM
    cdef list OE, T_ijk, r0, v0
    cdef double i, OMEGA, omega, mass
    cdef object date_initial

    OE = [initial_state.a, initial_state.e, initial_state.i, initial_state.OMEGA, initial_state.omega, initial_state.vi, initial_state.m]
    if (OE[0] > (3400 + 50
                 + 200) * 1e3) and (args.drag_passage == False) and args.body_shape == 'Spacecraft':
        index_steps_EOM = 3
    else:
        index_steps_EOM = 1
        # r_p = OE[0] * (1 - OE[1]) commented because not used 08-12
        #if args.res:
        #   print('Integration step number =', index_steps_EOM)

    # Rotation Matrix ijk tp pqw
    i = OE[2]
    OMEGA = OE[3]
    omega = OE[4]
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
    mass = OE[-1]

    # Clock
    date_initial = datetime(year=m.initialcondition.year,month=m.initialcondition.month,day=m.initialcondition.day,hour=m.initialcondition.hour,minute=m.initialcondition.min,second=m.initialcondition.second)

    config.count_numberofpassage = config.count_numberofpassage + 1
    # print('count_numberofpassage',config.count_numberofpassage)
    if config.count_numberofpassage != 1:
        t_prev = config.solution.orientation.time[-1]
    else:
        t_prev = m.initialcondition.time_rot
        # print('t_prev',t_prev)

    def f( t0,  in_cond,  m,  index_phase_aerobraking):
        ## Initialization
        cdef int passage_number, Mars_Gram_recalled_at_periapsis
        cdef object time_real, timereal
        cdef np.ndarray [double, ndim=1] pos_ii, vel_ii, pos_pp, vel_pp, pos_pp_hat, pos_ii_hat, vel_pp_hat, h_ii, h_pp, h_pp_hat
        cdef np.ndarray [double, ndim=1] uD, uE, uN, wind_pp, vel_pp_rw, vel_pp_rw_hat
        cdef double rot_angle
        cdef list omega_planet, OE, LatLong, uDuNuE, wind, y, L_PI
        cdef np.ndarray[double, ndim=1] gravity_ii, lift_pp_hat, drag_pp_hat, drag_pp, lift_pp, drag_ii, lift_ii,D_L_per_pp_hat
        cdef np.ndarray[double, ndim=1] thrust_pp_hat, thrust_pp, thrust_ii, force_ii
        cdef double mass, pos_ii_mag, vel_ii_mag,gamma,mu_fluid,area_tot, pos_pp_mag, vel_pp_mag, vi, h_ii_mag, h_pp_mag
        cdef double arg, gamma_ii, gamma_pp, lat, lon, alt, vN, vE, azi_pp, rho, T_p, p, length_car
        cdef double Re, sound_velocity, Mach, S, heat_load, Kn, heat_rate, cp, T_r
        cdef double aoa, wE, wN, wU, q
        cdef double bank_angle, CL, CD, beta, thrust_pp_mag, delta_v, energy

        ## Counters
        # Counter for all along the simulation of all passages
        config.count_aerobraking = config.count_aerobraking + 1
        passage_number = config.count_aerobraking
        # Counter for one entire passage
        config.count_dori = config.count_dori + 1
        # Counter for one phase
        config.count_phase = config.count_phase + 1


        # Clock
        time_real = date_initial + timedelta(seconds=t0)
        timereal = clock(time_real.year, time_real.month, time_real.day, time_real.hour, time_real.minute, time_real.second)

        # Assign states
        pos_ii = in_cond[0:3]  # Inertial position
        pos_ii += 0.
        vel_ii = in_cond[3:6]  # Inertial velocity
        vel_ii += 0.
        mass = in_cond[6]  # Mass, kg
        pos_ii_mag = np.linalg.norm(pos_ii)  # Inertial position magnitude
        vel_ii_mag = np.linalg.norm(vel_ii)  # Inertial velocity magnitude

        # Assign parameters
        omega_planet = m.planet.omega
        gamma = m.planet.gamma
        mu_fluid = m.planet.mu_fluid
        area_tot = m.body.Area_tot

        # TRANSFORM THE STATE
        # Inertial to planet relative transformation
        [pos_pp, vel_pp] = r_intor_p(pos_ii, vel_ii, m.planet, t0,t_prev) # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
        pos_pp_mag = np.linalg.norm(pos_pp)  # m, magnitude of position vector in PCPF
        pos_pp_hat = pos_pp / pos_pp_mag  # nd, unit position vector in PCPF
        pos_ii_hat = pos_ii / pos_ii_mag  # Inertial position vector direction

        vel_pp_mag = np.linalg.norm(vel_pp)
        vel_pp_hat = vel_pp / vel_pp_mag

        # Orbital Elements
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, m.planet)
        vi = OE[5]
        # print(math.degrees(vi))
        # if index_phase_aerobraking == 2:
        #     exit()
        Mars_Gram_recalled_at_periapsis = False
        if ( 0<vi<math.pi*0.5) and config.ascending_phase == False:
            config.ascending_phase = True
        elif (math.pi*0.5<=vi<= math.pi) and config.ascending_phase == True and args.body_shape == 'Blunted Cone':
            config.ascending_phase = False

        if config.ascending_phase == True and config.MarsGram_recall == False:
            config.atmospheric_data = {}

        # Angular Momentum Calculations
        h_ii = np.cross(pos_ii, vel_ii)  # Inertial angular momentum vector[m ^ 2 / s]
        cdef int index
        cdef double item
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
                h_pp[index] = 0.0
                index = index + 1
        h_pp_mag = np.linalg.norm(h_pp)
        h_pp_hat = h_pp / h_pp_mag


        # Inertial flight path angle
        arg = np.median([-1, 1, h_ii_mag / (pos_ii_mag * vel_ii_mag)])  # limit to[-1, 1]
        gamma_ii = math.acos(arg)
        if np.inner(pos_ii, vel_ii) < 0.0:
            gamma_ii = -gamma_ii


        # Relative flight path angle
        arg = np.median([-1, 1, h_pp_mag / (pos_pp_mag * vel_pp_mag)])  # limit to[-1, 1]
        gamma_pp = math.acos(arg)
        if np.inner(pos_pp, vel_pp) < 0.0:
            gamma_pp = -gamma_pp

        ## Derived Quantity Calculations


        # Compute latitude and longitude
        LatLong = rtolatlong(pos_pp, m.planet)
        lat = LatLong[1]
        lon = LatLong[2]
        alt = LatLong[0]
        # print(alt/10**3)
        if aerobraking_phase == 2 or aerobraking_phase == 0 :
            if (pos_ii_mag-m.planet.Rp_e-args.EI * 1e3) <= 0.0 and (config.drag_state == False) and config.ascending_phase == False:
                config.drag_state = True
                config.time_IEI = t0
            elif (pos_ii_mag-m.planet.Rp_e >= args.EI * 1e3) and (config.drag_state == True) and config.ascending_phase:
                config.drag_state = False
                config.time_OEI = t0

        if aerobraking_phase == 2 or aerobraking_phase == 0 :
            if args.control_mode == 1:
                x = 120
            else:
                x = 140

            if (config.heat_rate_prev>0.005 or abs(pos_ii_mag- m.planet.Rp_e <= x*1e3)) and (config.sensible_loads == False) and config.ascending_phase == False:
                config.sensible_loads = True
            elif (config.heat_rate_prev>0.005)  and (config.sensible_loads == True) and config.ascending_phase:
                config.sensible_loads = False

        # Compute NED basis unit vectors
        uDuNuE = latlongtoNED(LatLong)  # nd
        uD = uDuNuE[0]
        uE = uDuNuE[2]
        uN = uDuNuE[1]

        # compute azimuth
        vN = np.inner(vel_pp, uN)  # m / s
        vE = np.inner(vel_pp, uE)  # m / s
        azi_pp = math.atan2(vE, vN)  # rad

        # Get density, pressure, temperature and winds
        config.MarsGram_justrecalled=0
        rho,T_p,wind = dm(h=alt, p=m.planet, OE=OE, lat=lat, lon=lon, timereal=timereal, t0=t0, tf_prev = t_prev,
                          montecarlo=MonteCarlo,Wind=wind_m, args=args)

        # Define output.txt containing density data
        p = 0.0
        if args.body_shape == 'Spacecraft':
            length_car = m.body.length_SA + m.body.length_SC
        elif args.body_shape == 'Blunted Cone':
            length_car = m.body.BaseRadius *2
        Re = vel_pp_mag * rho * length_car / mu_fluid  # Reynolds number
        # Mach number
        sound_velocity = (gamma * m.planet.R * T_p) ** 0.5
        Mach = vel_pp_mag / sound_velocity
        S = ((gamma *0.5) ** 0.5) * Mach  # molecular speed ratio
        heat_load = in_cond[7]

        if config.drag_state == True:
            ## Check type of fluid
            Kn = 1.26 * (gamma) ** 0.5 * Mach / (Re+1e-5)
            if index_phase_aerobraking == 2:
                if (alt < 80000.0) & (config.index_warning_alt == 0):
                    if args.print_res:
                        print('WARNING: Altitude < 80 km!')
                    config.index_warning_alt = 1
                elif alt > 100000.0:
                    config.index_warning_alt = 0
                if (Kn < 0.1) and (config.index_warning_flow == 0):
                    if args.print_res:
                        print('Warning: transitional flow passage!')
                    config.index_warning_flow = 1
                elif Kn >= 0.1:
                    config.index_warning_flow = 0

        config.heat_load_past = heat_load
        # Heat rate and Control
        if (index_phase_aerobraking == 2 or index_phase_aerobraking == 1.75 or index_phase_aerobraking == 2.25) and config.drag_state and config.initial_position_closed_form:
            #evaluates the closed form solution the first time at EI km
            if abs(pos_ii_mag-m.planet.Rp_e-args.EI*1e3) <=1e-2 and (args.control_mode==2 or args.control_mode==3) and config.time_switch_1 == 0:
                cntr(ip, m , t=(t0 - config.time_IEI) , position=config.initial_position_closed_form ,
                                  args=args ,
                                  index_ratio=[1 , 0] , state=[T_p , rho , S])

            if config.MarsGram_justrecalled == True and config.index_Mars_Gram_call != 1:  # in MC, when we reavaluate Mars Gram there is a discontinuity with the density which is created by how the density data are created. This discontinuity create really high peaks. We recalculate aoa for the new density data.
                config.aoa = cntr(ip, m , t=(t0 - config.time_IEI) , position=config.initial_position_closed_form ,
                                  args=args ,
                                  index_ratio=[1 , 0] , state=[T_p , rho , S])
                config.state_flesh1.append([T_p, rho, S])
                heat_rate = hr_c(S, T_p, m, rho, vel_pp_mag,config.aoa)


            if (index_phase_aerobraking == 2):# and args.control_mode == 3) or (args.control_mode == 2):
                if args.control_in_loop:  # Control in-loop online
                    config.state_flesh1 = [T_p,rho,S]
                    config.aoa = cntr(ip, m , t=(t0 - config.time_IEI) , position=config.initial_position_closed_form , args=args ,
                               index_ratio=[1 , 1] , state=config.state_flesh1 , current_position=OE)
                elif args.control_in_loop == False and args.integrator == 'Python':
                    if config.ctrller.count_controller != config.ctrller.count_prev_controller and config.ctrller.stored_state==0 and t0!=config.ctrller.prev_time:
                        # print('t integrator',t0)
                        config.state_flesh1.append([T_p , rho , S])
                        if config.ctrller.count_controller == 2:
                            state = config.state_flesh1[-1]
                        else:
                            state = config.state_flesh1[-2]
                            config.state_flesh1.pop(0)
                        config.ctrller.stored_state = 1
                        config.ctrller.prev_time = time_0
                        config.aoa = cntr(ip, m , t=(t0-config.time_IEI) , position=config.initial_position_closed_form , args=args ,
                              index_ratio=[1,1], state= state, current_position = OE)
            # Heat Rate
            heat_rate = hr_c(S, T_p, m, rho, vel_pp_mag,config.aoa)
            cp = m.planet.gamma / (m.planet.gamma - 1.0) * m.planet.R

            T_r = 0.0
            #emissivity = 0.82
            #stef_bol_const = 5.67*10**-4
            #config.T_w = (heat_rate / (emissivity*stef_bol_const))**-4 - 4

        else:
            T_r = 0.0
            heat_rate = 0.0

        aoa = config.aoa
        config.heat_rate_prev = heat_rate # save current heat rate

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
        rot_angle = np.linalg.norm(omega_planet) * (t0+t_prev)  # rad
        L_PI = [[math.cos(rot_angle), math.sin(rot_angle), 0.0],
                [-math.sin(rot_angle), math.cos(rot_angle), 0.0],
                [0.0, 0.0, 1.0]]

        L_PI = [[x+0. for x in row] for row in L_PI]

        gravity_ii = mass * gm(pos_ii_mag, pos_ii, p=m.planet, mass=mass, vel_ii=vel_ii)

        # Vectors and Rotation
        # Tensors of interest


        # n1 = np.array([[0, h_pp_hat[2], - h_pp_hat[1]],
        #                [-h_pp_hat[2], 0, h_pp_hat[0]],
        #                [h_pp_hat[1], -h_pp_hat[0], 0]])
        #
        # R1 = np.eye(3) + math.sin(math.pi/2) * n1 + (1 - math.cos(math.pi/2)) *n1*n1#* np.inner(n1, n1)
        #
        # #print(n1*n1)
        # #print(np.inner(n1, n1))
        # n2 = np.array([[0, -vel_pp_rw_hat[2],  vel_pp_rw_hat[1]],
        #                [vel_pp_rw_hat[2], 0, -vel_pp_rw_hat[0]],
        #                [-vel_pp_rw_hat[1], vel_pp_rw_hat[0], 0]])  # change with vel_pp_rw_hat
        # bank_angle = 0
        # R2 = np.eye(3) + math.sin(bank_angle) * n2 + (1 - math.cos(bank_angle)) *n2*n2# np.inner(n2, n2)


        bank_angle = 0.0
        lift_pp_hat = np.cross(h_pp_hat,vel_pp_rw_hat) #perpendicular vector to angular vector and velocity

        # Vehicle Aerodynamic Forces
        # CL and CD
        [CL,CD] = am(aoa=aoa, body=m.body, T=T_p, S=S, args=args, montecarlo=MonteCarlo)
        beta = mass/(CD*area_tot) # ballistic coefficient
        # Force calculations
        drag_pp_hat = -vel_pp_rw_hat # Planet relative drag force direction

        drag_pp = q * CD * area_tot * drag_pp_hat  # Planet relative drag force vector
        lift_pp = q * CL * area_tot * lift_pp_hat* math.cos(bank_angle)  # Planet relative lift force vector

        drag_ii = np.inner(np.transpose(L_PI), drag_pp)  # Inertial drag force vector
        lift_ii = np.inner(np.transpose(L_PI), lift_pp)  # Inertial lift force vector

        # Check if propellant mass is greater than 0 kg
        if config.index_propellant_mass == 1:
            if mass-args.dry_mass <= 0.5:
                config.index_propellant_mass = 0
                m.engine.T = 0
                if args.print_res:
                    print('No Fuel Left!')

        # Odyssey_firing_plan(delta_v , args)

        # Thrust
        delta_v = (m.engine.g_e * m.engine.Isp) * np.log(initial_state.m / mass)
        thrust_pp_mag = thrust_m(t0,m.engine.T,delta_v,args,index_phase_aerobraking) # N

        # Rodrigues rotation formula to rotate thrust vector of angle phi around angular vector from D direction
        D_L_per_pp_hat = (np.cross(drag_pp_hat,lift_pp_hat))
        thrust_pp_hat = drag_pp_hat*math.cos(args.phi) + np.cross(D_L_per_pp_hat,drag_pp_hat) *math.sin(args.phi)+D_L_per_pp_hat*(np.inner(D_L_per_pp_hat,drag_pp_hat))*(1-math.cos(args.phi))
        #these two ways give the same direction
        #thrust_pp_hat = -vel_pp_hat*math.cos(args.phi)+np.cross(h_pp_hat,-vel_pp_hat)*math.sin(args.phi)+h_pp_hat*(np.inner(h_pp_hat,-vel_pp_hat))*(1-math.cos(args.phi))
        # thrust_pp_hat = vel_pp_hat*math.cos(args.phi)
        thrust_pp = thrust_pp_mag * thrust_pp_hat#; % N
        thrust_ii = np.inner(np.transpose(L_PI), thrust_pp)

        # Total Force
        # Total inertial external force vector on body [N]
        force_ii = drag_ii + lift_ii + gravity_ii + thrust_ii

        # EOM
        ydot = np.zeros(8)

        ydot[0:3] = vel_ii
        ydot[3:6] = force_ii / mass
        ydot[6] = -np.linalg.norm(thrust_ii) / (m.engine.g_e * m.engine.Isp) # mass variation
        ydot[7] = heat_rate
        energy = (vel_ii_mag ** 2)*0.5- m.planet.mu / (pos_ii_mag)

        # ## SAVE RESULTS
        if config.results_save:
            sol = np.hstack([[t0], [timereal.year], [timereal.month], [timereal.day], [timereal.hour],[timereal.minute],
                         [timereal.second], [numberofpassage], pos_ii, vel_ii, [pos_ii_mag], [vel_ii_mag], pos_pp,
                         [pos_pp_mag], vel_pp, [vel_pp_mag], [OE[0]], [OE[1]], [OE[2]], [OE[3]], [OE[4]], [OE[5]],
                         [lat], [lon], [alt], [gamma_ii], [gamma_pp], h_ii, h_pp, [h_ii_mag], [h_pp_mag], uD, uE, uN,
                         [vN], [vE], [azi_pp], [rho], [T_p], [p], wind, [CL], [CD], [aoa], [S], [mass], [heat_rate], [heat_load],[T_r], [q],
                         gravity_ii, drag_pp, drag_ii, lift_pp, lift_ii, force_ii, [energy], config.index_MonteCarlo, int(config.drag_state)]).reshape((1, 91)).tolist()
            config.solution_intermediate.extend(sol)
        return ydot

    ## EVENTS

    def eventfirststep(double t, object y):
        # print(((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - 250*1e3))
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - 250*1e3)
    eventfirststep.direction = -1
    eventfirststep.terminal = True

    def eventsecondstep(double t, object y):
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - 260*1e3)

    eventsecondstep.terminal = True
    eventsecondstep.direction = 1

    def reached_EI(double t, object y):
        return (y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - args.EI * 1e3
    reached_EI.terminal = False
    reached_EI.direction = -1

    def reached_AE(double t, object y):
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - args.AE*1e3)
    reached_AE.terminal = False
    reached_AE.direction = 1

    def out_drag_passage(double t, object y):
        if abs(((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - args.AE*1e3))  <= 1e-5:
            if args.heat_load_sol == 0 or args.heat_load_sol == 2:
                config.aoa = m.aerodynamics.aoa
            elif args.heat_load_sol == 1 or args.heat_load_sol == 3:
                config.aoa = 0.0
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - args.AE*1e3)

    out_drag_passage.terminal = True
    out_drag_passage.direction = 1

    def in_drag_passage(double t, object y):
        cdef list pos_ii, vel_ii, OE_closedform
        pos_ii = [y[0], y[1], y[2]]  # Inertial position
        vel_ii = [y[3], y[4], y[5]]  # Inertial velocity
        [pos_pp, vel_pp] = r_intor_p(pos_ii, vel_ii, m.planet, t,
                                     t_prev)  # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
        LatLong = rtolatlong(pos_pp, m.planet)
        h0 = LatLong[0]
        # Control Setting
        if ip.gm == 2:
            cond = h0- args.EI*1e3
            thr = 500
        else:
            cond = (y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - args.EI*1e3
            thr = 1e-5
        if abs(cond) <= thr:
        # if abs((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - args.EI*1e3) <= 1e-5:
            config.ctrller.guidance_t_eval = np.arange(t , t+1500 , 1 / args.flash1_rate).tolist()

            # State definition for control 2,3 State used by the closed-form solution
            pos_ii = [y[0] , y[1] , y[2]]  # Inertial position
            vel_ii = [y[3] , y[4] , y[5]]  # Inertial velocity
            OE_closedform = rvtoorbitalelement(pos_ii , vel_ii , y[6] , m.planet)
            config.initial_position_closed_form = OE_closedform
        return cond
    in_drag_passage.terminal = True
    in_drag_passage.direction = -1

    def in_drag_passage_nt(double t, object y):
        cdef list pos_ii, vel_ii, OE_closedform
        # Control Setting
        # print(config.initial_position_closed_form)
        pos_ii = [y[0], y[1], y[2]]  # Inertial position
        vel_ii = [y[3], y[4], y[5]]  # Inertial velocity
        [pos_pp, vel_pp] = r_intor_p(pos_ii, vel_ii, m.planet, t,
                                     t_prev)  # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
        LatLong = rtolatlong(pos_pp, m.planet)
        h0 = LatLong[0]
        if ip.gm == 2 and not args.drag_passage:
            cond = (h0- args.EI*1e3)
            thr = 500
        else:
            cond = ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - args.EI*1e3)
            thr = 1e-3
        if abs(cond) <= thr and not config.initial_position_closed_form:
            config.ctrller.guidance_t_eval = np.arange(t , t+1500 , 1 / args.flash1_rate).tolist()

            # State definition for control 2,3 State used by the closed-form solution
            OE_closedform = rvtoorbitalelement(pos_ii , vel_ii , y[6] , m.planet)
            config.initial_position_closed_form = OE_closedform
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - args.EI*1e3)
    in_drag_passage_nt.terminal = False
    in_drag_passage_nt.direction = -1


    def apoasispoint(double t, object y):
        cdef double vi
        cdef list pos_ii, vel_ii
        pos_ii = [y[0] , y[1] , y[2]]  # Inertial position
        vel_ii = [y[3] , y[4] , y[5]]  # Inertial velocity
        vi = rvtoorbitalelement(pos_ii , vel_ii , y[6] , m.planet)[5]
        return math.degrees(vi)-180.0

    apoasispoint.terminal = True
    apoasispoint.direction = 1

    def periapsispoint(double t, object y):
        cdef double vi
        cdef list pos_ii, vel_ii
        pos_ii = [y[0] , y[1] , y[2]]  # Inertial position
        vel_ii = [y[3] , y[4] , y[5]]  # Inertial velocity
        vi = rvtoorbitalelement(pos_ii , vel_ii , y[6] , m.planet)[5]
        return vi

    periapsispoint.terminal = False
    periapsispoint.direction = 1

    def impact(double t, object y):  # nota: Event impact has to be the last because of the impact condition in the aerobraking script.
        cdef double min_alt
        if args.body_shape == 'Blunted Body':
            min_alt = 1.0*10**3
        else:
            min_alt = 35.0*10**3
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - (
                    m.planet.Rp_e + min_alt))  # consider impact when the altitude is less than 35 km

    impact.terminal = True

    def apoapsisgreaterperiapsis(double t, object y):
        cdef object r,v
        cdef np.ndarray[double, ndim=1] h
        cdef double Energy, e
        cdef double a, r_a, r_p
        y = [x+0. for x in y]
        r = y[0:3]
        v = y[3:6]
        Energy = (np.linalg.norm(v)**2)*0.5 - m.planet.mu / (np.linalg.norm(r))
        a = - m.planet.mu / (2 * Energy)
        h = np.cross(r, v)
        h += 0.
        e = (1 + (2.0 * Energy * (np.inner(h, h)) / (m.planet.mu) ** 2)) ** 0.5
        # print('e2',np.linalg.norm(np.multiply(1/m.planet.mu,np.multiply((np.linalg.norm(v)**2-m.planet.mu/np.linalg.norm(r)),r)-np.multiply((np.dot(v,r)),v))))

        r_a = a * (1.0 + e)
        r_p = a * (1.0 - e)
        if r_a < (r_p) and args.body_shape == 'Spacecraft':
            print('Periapsis greater than apoapsis!')
        return (r_a - r_p)
    apoapsisgreaterperiapsis.terminal = True

    def stop_firing(double t, object y):
        cdef double mass, delta_v
        mass = y[6]
        delta_v = (m.engine.g_e * m.engine.Isp)*np.log(initial_state.m/mass) # delta-v variation
        return delta_v - args.delta_v#(delta_v_max+config.delta_v_man)
    stop_firing.terminal = True

    def guidance(double t, object y):
        # print(t - (t_prev + 1 / args.trajectory_rate))
        if config.ctrller.stored_state == 1:
            config.ctrller.count_prev_controller = config.ctrller.count_controller
        if t - config.ctrller.t >1:
            print('Decrease step size of integration and tolerance')
        if abs(t-(config.ctrller.t)) <= 1e-08: # if give you error decreases tolerance and decreases this value
            config.ctrller.stored_state = 0
            config.ctrller.count_controller += 1
            # print(config.ctrller.count_controller,config.ctrller.guidance_t_eval[config.ctrller.count_controller])
            config.ctrller.t = config.ctrller.guidance_t_eval[config.ctrller.count_controller]
        return t-(config.ctrller.t)
    guidance.terminal = False
    guidance.direction = 1

    def heat_rate_check(double t, object y):
        cdef int x
        if args.control_mode == 1:
            x = 120
        else:
            x = 140
        if abs((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - x*1e3)<= 1e-5:
            config.ctrller.guidance_t_eval = np.arange(t , t+2500 , 1 / args.flash1_rate).tolist()
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - x*1e3)
    heat_rate_check.terminal = True

    def heat_load_check_exit(double t, object y):
        cdef int x
        if args.control_mode == 1:
            x = 120
        else:
            x = 160
        return ((y[0] ** 2 + y[1] ** 2 + y[2] ** 2) ** 0.5 - m.planet.Rp_e - x * 1e3)
    heat_load_check_exit.terminal = True
    heat_load_check_exit.direction = 1

    cdef int stop_simulation, save_pre_index, save_post_index, range_phase_i, aerobraking_phase, i_sim, final_conditions_notmet, min_index
    cdef list in_cond, events, time_solution, y
    cdef float index_phase_aerobraking, save_ratio
    cdef double step, time_0
    cdef str simulator, method
    cdef double r_tol, a_tol, length_sim, initial_time, final_time, t, min_value
    cdef object solution
    cdef list pos_ii_gnc, vel_ii_gnc, OE_gnc, state_flesh1

    time_0 = 0
    if args.heat_load_sol == 0 or args.heat_load_sol == 2:
        config.aoa = m.aerodynamics.aoa
    elif args.heat_load_sol == 1 or args.heat_load_sol == 3:
        config.aoa = 0.0
    stop_simulation = False
    save_pre_index = 0
    save_post_index = 0
    config.impact = False
    config.solution_intermediate = []
    config.count_dori = 0
    config.atmospheric_data = {}
    config.previous_atmospheric_data = {}
    config.ascending_phase = False
    config.evaluate_switch_heat_load = False
    config.state_inner_boundary_atmosphere = [] # used in Density model for vi def
    config.time_IP = config.time_OP
    config.time_IEI = 0
    config.time_OEI = 0
    config.time_switch_1 = 0.0
    config.time_switch_2 = 0.0
    if args.heat_load_sol == 2:
        config.time_switch_2 = 1000.0
    config.timer_revaluation = 0
    config.closed_form_solution_off = 1 # used in closed form solution online to run the solution only once
    config.initial_position_closed_form = []
    config.heat_rate_list = []  # checked this - if used
    config.aoa_list = []  # checked this - if used
    config.ctrller.guidance_t_eval = []
    config.ctrller.count_prev_controller = 0
    config.ctrller.count_controller = 1
    config.ctrller.stored_state = 1
    config.ctrller.prev_time = 0
    config.ctrller.t = 0
    config.security_mode = False
    config.stop_simulation=False
    config.results_save = 1
    config.drag_state = False
    config.aoa_past = m.aerodynamics.aoa


    if (r0[0]**2+r0[1]**2+r0[2]**2)**0.5-m.planet.Rp_e <= args.EI*1e3:
        config.drag_state = True
        config.initial_position_closed_form = OE

    config.sensible_loads = False

    config.counter_integrator =0


    # Def initial conditions
    in_cond = [r0[0], r0[1], r0[2], v0[0], v0[1], v0[2],mass,0]

    # If aerobraking maneuver allowed, add a prephase 0
    range_phase_i = 1
    if args.thrust_control == 'Aerobraking Maneuver':
        range_phase_i = 0

    # Solve Equations of Motion
    for aerobraking_phase in range(range_phase_i,4):
        index_phase_aerobraking = aerobraking_phase
        if (index_steps_EOM == 1 or args.drag_passage) and (aerobraking_phase == 1 or aerobraking_phase == 3 or aerobraking_phase == 0):
            continue

    # Definition of events
        if aerobraking_phase == 0:
            events = [stop_firing , apoapsisgreaterperiapsis , impact]
        elif aerobraking_phase == 1:
            events = [eventfirststep, apoapsisgreaterperiapsis , impact]
        elif aerobraking_phase == 2 and args.drag_passage:  # this conditions is when control mode !=0
            events = [out_drag_passage, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis , impact]
        elif aerobraking_phase == 2 and index_steps_EOM == 1 and args.body_shape == 'Blunted Cone':
            events = [out_drag_passage, apoasispoint , periapsispoint , impact]
        elif aerobraking_phase == 2 and index_steps_EOM == 1 and args.drag_passage == False:
            events = [apoasispoint, periapsispoint, in_drag_passage_nt, apoapsisgreaterperiapsis , impact]
        elif aerobraking_phase == 2:  # this condition is when control mode == 0
            events = [eventsecondstep, periapsispoint, reached_EI,reached_AE, in_drag_passage_nt, apoapsisgreaterperiapsis, impact]
        elif aerobraking_phase == 3:
            events = [apoasispoint , periapsispoint , apoapsisgreaterperiapsis, impact]

        if index_phase_aerobraking == 2:
            save_pre_index = len(config.solution.orientation.time)
            simulator = args.integrator


    # Definition step size
        if aerobraking_phase == 1 or aerobraking_phase == 3:
            step = 5#10 # real value 5
            r_tol = 1e-10#1e-4 # real value 1e-10
            a_tol = 1e-12#1e-5 # real value 1e-12
            simulator = 'Python'
            method = 'BDF'#'RK23'#real value 'BDF'
            save_ratio = 5
        elif aerobraking_phase == 0:
            step = 0.1
            r_tol = 1e-9
            a_tol = 1e-11
            simulator = 'Python'
            method = 'BDF'
            save_ratio = 5
        elif aerobraking_phase == 2:
            if args.integrator == 'Python':
                step = 0.1#1 # real value 0.1 Change back to the commented values
                r_tol = 1e-9#1e-5 # real value 1e-9
                a_tol = 1e-11#1e-7 # real value 1e-11
                if MonteCarlo:
                    method = 'RK45'
                else:
                    method = 'BDF'
                save_ratio = 1
            else:
                step = 1/(args.trajectory_rate)
                save_ratio = round(args.trajectory_rate / args.save_rate)
        # Definition length simulation:
        length_sim = 1e10
        i_sim = 0
        time_solution = []

        config.continue_simulation = True


        while config.continue_simulation:
            index_phase_aerobraking = aerobraking_phase
            # if control mode =! 0, redefine sim setting and creates two more phases until reaching EI and out of the AE
            #phase 2: between 120 km alt
            if aerobraking_phase == 2 and (args.control_mode != 0 and args.control_in_loop == 0 and config.drag_state == True and config.sensible_loads == True and config.ascending_phase == False):
                simulator = args.integrator
                events = [out_drag_passage, heat_load_check_exit, periapsispoint , apoapsisgreaterperiapsis, impact]
                #count += 1
                if args.integrator == 'Costumed':
                    step = 1/(args.trajectory_rate)
                else:
                    method = 'RK23'
                    if args.control_in_loop: # online control
                        step = 10/(args.trajectory_rate)#0.5
                    else: # more realistic control: use previous state and control rate != trajectory rate
                        step = 1/(4*args.flash1_rate) # don't change to higher values otherwise controller doesn't work0.05#
                        r_tol = 1e-6
                        a_tol = 1e-7
                        if step>= args.flash1_rate:
                            step = 1/(2*args.flash1_rate)
                        events.insert( 3, guidance)
                config.ctrller.t = config.ctrller.guidance_t_eval[config.ctrller.count_controller]
            #phase 1.75: between EI km alt and 120 km
            elif aerobraking_phase == 2 and args.control_mode != 0 and args.control_in_loop == 0 and config.ascending_phase == False and config.drag_state == True:  # prephase 1 between EI km and heat loads control
                events = [periapsispoint , out_drag_passage, heat_rate_check,apoapsisgreaterperiapsis , impact]
                simulator = 'Python'
                index_phase_aerobraking = 1.75
                step = 0.05#0.5
                r_tol = 1e-6#1e-7
                a_tol = 1e-7#1e-8
                method = 'RK23'
            # phase 1.5: between 250 km alt and EI km
            elif aerobraking_phase == 2 and args.control_mode != 0 and args.control_in_loop == 0  and config.ascending_phase == False and config.drag_state == False: # prephase 1 between 250 km to EI
                events = [periapsispoint, in_drag_passage,apoapsisgreaterperiapsis,impact]  # the terminal event is in_drag_passage
                simulator = 'Python'
                step = 0.1
                r_tol = 1e-9
                a_tol = 1e-10
                index_phase_aerobraking = 1.5
                method = 'RK23'
            # phase 2.25: between 120 km alt and AE km
            elif aerobraking_phase == 2 and args.control_mode != 0 and args.control_in_loop == 0  and config.ascending_phase and config.drag_state == True:
                events = [periapsispoint, out_drag_passage,apoapsisgreaterperiapsis,impact]
                simulator = 'Python'
                index_phase_aerobraking = 2.25
                step = 0.1
                r_tol = 1e-9
                a_tol = 1e-10
                method = 'RK23'
            # phase 2.5: between AE km alt and 250 km
            elif aerobraking_phase == 2 and args.control_mode != 0 and args.control_in_loop == 0 and config.ascending_phase:
                events = [eventsecondstep , periapsispoint, eventsecondstep,apoapsisgreaterperiapsis,impact]
                simulator = 'Python'
                index_phase_aerobraking = 2.5
                step = 0.5
                r_tol = 1e-8
                a_tol = 1e-9
                method = 'RK23'

            if args.print_res:
                print('Step #' , index_phase_aerobraking)

            if simulator == 'Python':
                ## Python Integrator
                # Time initialization
                initial_time , final_time = time_0, time_0 + length_sim

                # Run Simulation
                eomfunction = lambda t , in_cond: f(t , in_cond , m , index_phase_aerobraking)
                solution = solve_ivp(eomfunction , t_span=(initial_time , final_time) , y0=in_cond ,max_step=step, method=method ,
                                 dense_output=True  , rtol=r_tol,atol=a_tol, events=events)  # (DOP853) , rtol=1e-9, atol=1e-11
                config.counter_integrator += 1
                in_cond = [solution.y[0][-1] , solution.y[1][-1] , solution.y[2][-1] , solution.y[3][-1] ,
                           solution.y[4][-1] ,
                           solution.y[5][-1] , solution.y[6][-1] , solution.y[7][-1]]
                # Save results
                time_solution.extend(solution.t)
                time_0 = time_solution[-1]
                if aerobraking_phase == 0:
                    from utils.misc import new_periapsis
                    new_periapsis(m, in_cond[0:3], in_cond[3:6], args)

            elif simulator == 'Costumed':
                if args.integrator == 'Costumed':
                    class solution():
                        def __init__ ( self ):
                            self.t_events = [[] , []]

                    solution = solution()

                while stop_simulation == False:
                    initial_time = time_0
                    ## Costumed Integrator
                    config.MarsGram_recall = 1  # MarsGram can be recalled only when a new step calculation just beginned. This is to avoid to have inconsistency  in the atmospheric data for the same time step
                    config.results_save = 0
                    y , t , stop_simulation, solution = RK4(f , step , initial_time , in_cond , m , T_ijk ,
                                                                  index_phase_aerobraking , args, solution)

                    in_cond = [y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7]]
                    config.counter_integrator += 1
                    config.results_save = 1
                    f(t , np.array(in_cond) , m , index_phase_aerobraking)
                    # New initial condition

                    time_0 = t

                    # Save results
                    time_solution.append(t)

                    ## Guidance, Navigation and Control
                    if args.control_mode != 0:
                        if (i_sim % (round(args.trajectory_rate/args.flash1_rate)) == 0) and not args.control_in_loop:
                            temp = y[0:3]  # Inertial position
                            pos_ii_gnc = [value+0. for value in temp]
                            temp = y[3:6]  # Inertial velocity
                            vel_ii_gnc = [value+0. for value in temp]
                            OE_gnc = rvtoorbitalelement(pos_ii_gnc , vel_ii_gnc , y[6] , m.planet)
                            state_flesh1 = [config.solution_intermediate[-1][56] ,config.solution_intermediate[-1][55] , config.solution_intermediate[-1][64]]
                            config.aoa = cntr(ip, m , t=(t - config.time_IEI) , current_position=OE ,
                                          position=config.initial_position_closed_form , args=args ,
                                          index_ratio=[1,1] , state=state_flesh1)

                i_sim += 1

            # Define breaker campaign impact or apoapsis greater than periapsis
            continue_campaign = event(solution)
            if continue_campaign == False:
                config.impact = True
                break

            # Breaker conditions
            if simulator == 'Python':
                if len(solution.t_events[0]) != 0 or (args.drag_passage and index_phase_aerobraking == 2.25 and len(solution.t_events[1]) != 0):
                    config.continue_simulation = False
                    break
            else:
                if args.control_mode == 0 and stop_simulation == True:
                    config.continue_simulation = False
                    break

        # Save Results
        #round(args.trajectory_rate / args.save_rate)
        time_0 = save_results(time_solution , save_ratio)

        if index_phase_aerobraking == 2 or index_phase_aerobraking == 2.5 or (index_phase_aerobraking == 2.25 and args.drag_passage):
            save_post_index = len(config.solution.orientation.time)


        # Re-Set count index to 0
        config.count_phase = 0

        # Define breaker campaign
        if continue_campaign == False:
            if save_post_index == 0:
                save_post_index = len(config.solution.orientation.time)
            break

    # Re-Set count index to 0
    config.count_dori = 0
    if args.drag_passage == False and (math.pi -
                config.solution.orientation.oe[-1][-1] ) > 1e-4 and continue_campaign == True and args.body_shape != 'Blunted Cone':
        final_conditions_notmet = True
        events = [apoasispoint]
    else:
        final_conditions_notmet = False
    count_temp = 0
    while final_conditions_notmet:
        in_cond = [config.solution.orientation.pos_ii[0][-1] , config.solution.orientation.pos_ii[1][-1] , config.solution.orientation.pos_ii[2][-1] ,
                       config.solution.orientation.vel_ii[0][-1] ,config.solution.orientation.vel_ii[1][-1] ,config.solution.orientation.vel_ii[2][-1] ,
                       config.solution.performance.mass[-1] , config.solution.performance.heat_load[-1]]
        initial_time , final_time = time_0, time_0 + 100
        step = 0.005#0.5# real value 0.005
        r_tol = 1e-12#1e-7# real value 1e-12
        a_tol = 1e-13#1e-8# real value 1e-13

        # Run Simulation
        try:
            eomfunction = lambda t , in_cond: f(t , in_cond , m , index_phase_aerobraking)
            solution = solve_ivp(eomfunction , t_span=(initial_time , final_time) , y0=in_cond , max_step=step ,
                                 method='BDF' ,
                                 dense_output=True , rtol=r_tol , atol=a_tol ,
                                 events=events)  # (DOP853) , rtol=1e-9, atol=1e-11
            config.counter_integrator += 1
            time_0 = save_results(solution.t , ratio = 0.1)
            count_temp += 1
        except: #Exception:
        #     print(math.pi -config.solution.orientation.oe[-1][-1] )
        #     traceback.print_exc()
        #     print('Braked Final Conditions Loop')
            break
        if count_temp >5:
            break

        if args.drag_passage == False and (math.pi -
                config.solution.orientation.oe[-1][-1] ) > 1e-4 and continue_campaign == True:
            final_conditions_notmet = True
            events = [apoasispoint]
        else:
            final_conditions_notmet = False

    config.save_index_heat = len(config.solution.orientation.time)
    config.time_OP = len(config.solution.orientation.time)

    config.altitudeperiapsis.append(min(config.solution.orientation.alt[save_pre_index:save_post_index])*1e-3)
    config.max_heatrate.append(max(config.solution.performance.heat_rate[save_pre_index:save_post_index]))
    config.delta_v_man = ((m.engine.g_e * m.engine.Isp)*np.log(m.body.Mass/config.solution.performance.mass[-1]))

    if args.print_res:
        try:
            #print("Thermal limit overcomed {} times!".format(count))
            print("Actual periapsis altitude {0:.2f}km - Vacuum periapsis altitude = {1:.2f}km".format(min(config.solution.orientation.alt[save_pre_index:save_post_index])*1e-3, (config.solution.orientation.oe[0][-1] * (1 - config.solution.orientation.oe[1][-1]) - m.planet.Rp_e)*1e-3))
            print("Ra new = {0:.2f}km".format((config.solution.orientation.oe[0][-1] * (1 + config.solution.orientation.oe[1][-1]))*1e-3))
            print("Ra new = {}km".format(
                (config.solution.orientation.oe[0][-1] * (1 + config.solution.orientation.oe[1][-1])) * 1e-3))
            print('HEAT RATE IS',np.max(config.solution.performance.heat_rate[save_pre_index:save_post_index-1]),'W/cm^2')
            print('HEAT LOAD IS',config.solution.performance.heat_load[save_post_index-1],'J/cm^2')
            print('Fuel Mass is',(config.solution.performance.mass[-1]-args.dry_mass),' kg')
            print('Total time is', config.solution.orientation.time[save_post_index-1]-config.solution.orientation.time[save_pre_index],'s')
            print('Delta-v is',config.delta_v_man ,' m/s' )
            print('Delta-E is',(config.solution.forces.energy[config.time_OP-1]-config.solution.forces.energy[config.time_IP])*1e-3,'kJ')
            min_value = min(config.solution.orientation.alt[save_pre_index:save_post_index])
            min_index = config.solution.orientation.alt.index(min_value)
            print('Latitude of periapsis {}deg'.format(math.degrees(config.solution.orientation.lat[min_index])))
            print('Longitutde of periapsis {}deg'.format(math.degrees(config.solution.orientation.lon[min_index])))
            if args.body_shape == 'Blunted Cone':
                max_value = max(config.solution.performance.q)
                max_index = config.solution.performance.q.index(max_value)
                print('Max Dynamic Pressure {}N/m^2 at time {}s'.format(max_value,config.solution.orientation.time[max_index]))
                max_value = max(config.solution.performance.heat_rate)
                max_index = config.solution.performance.heat_rate.index(max_value)
                print('Max Heat Rate {}W/cm^2 at time {}s'.format(max_value ,
                                                                    config.solution.orientation.time[max_index]))
                max_value = max(config.solution.performance.heat_load)
                max_index = config.solution.performance.heat_load.index(max_value)
                print('Max Heat Load {}J/cm^2 at time {}s'.format(max_value , config.solution.orientation.time[max_index]))
        except:
            print('Problem in the indexes')
    # Results to delete, only for plot for odyssey mission
    config.periapsis_list.append(min(config.solution.orientation.alt[save_pre_index:save_post_index])*1e-3) # actual periapsis
    config.orbit_number_list.append(config.count_numberofpassage+1)
    config.delta_v_list.append(config.delta_v_man)

    return continue_campaign
