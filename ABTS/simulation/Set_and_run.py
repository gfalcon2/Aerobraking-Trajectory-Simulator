#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

# Main file to define purposes, principal loop, input files, graph and results
from physical_models.Mission import *
from simulation.Aerobraking import aerobraking
import config as cnf
import math
from physical_models.Planet_data import planet_data
import numpy as np
import time as timecode
import random
from datetime import *
import os
from utils.save_cvs import *
import time

def aeroraking_campaign(args,state):


    save_res = args.results


    # Descent towards Mars
    def mission_conf ():
        mission = {'Purpose': purpose , 'Gravity Model': args.gravity_model , 'Density Model': args.density_model , 'Wind': args.wind,
                   'Aerodynamic Model': args.aerodynamic_model , 'Thermal Model': args.thermal_model, 'Control': args.control_mode, 'Firings': args.thrust_control, 'Shape':args.body_shape,'Monte Carlo': args.montecarlo}
        return (mission)


    # First phase
    purpose = 'Aerobraking around Mars.'


    mission = mission_conf()
    if args.print_res == True:
        ('mission is', mission)

    ip = missiondef(mission_conf())
    p_class = planet_data(ip.M.planet)
    if args.gravity_model == 'Inverse Squared':
        p_class.Rp_p = p_class.Rp_e

    # Vehicle - calculation notebook page1
    # Mass
    dry_mass = args.dry_mass
    prop_mass = args.prop_mass
    mass = dry_mass + prop_mass

    # Spacecraft Shape
    if args.body_shape == 'Spacecraft':
        area_body = 7.26#33.38#7.26# This is recalculated for the new sc config. 11 (look notes)# m^2 2001 Mars Odyssey Aerobraking, Smith & Bell paper
        length_sp = 3.7617#11.4#3.7617#5.7 # m solar array length https://www.jpl.nasa.gov/news/press_kits/odysseyarrival.pdf
        height_sp = area_body/length_sp
        # Main Body
        length_ody = 2.2 # m
        height_ody = 1.7 # m
        width_ody =  2.6 # m
    #Blunted Body Shape
    elif args.body_shape == 'Blunted Cone':
        delta = 70#45 #deg
        nose_radius = 0.6638#0.3 #m
        base_radius = 2.65/2#0.5 #m


    apoapsis = state['Apoapsis']

    state['Periapsis'] = p_class.Rp_e + (state['Periapsis'])*1e3 #98
    state['vi'] = np.radians(180.0001)

    if args.montecarlo == True:
        from physical_models.MonteCarlo_perturbations import monte_carlo_initial_conditions
        state = monte_carlo_initial_conditions(state,args)



    # Initial Condition Calcs
    semimajoraxis_in = (state['Apoapsis'] + state['Periapsis'])/2
    eccentricity_in = (state['Apoapsis'] - state['Periapsis'])/ (state['Apoapsis'] + state['Periapsis'])
    apoapsis = state['Apoapsis']
    periapsis = state['Periapsis']
    # Initial Condition
    if args.drag_passage == True:
        h_0 = 160 * 10 ** 3
    elif args.body_shape == 'Blunted Cone':
        h_0 = 120*10**3
        args.AE = args.EI = h_0/10**3
    if args.drag_passage or args.body_shape == 'Blunted Cone':
        r = p_class.Rp_e + h_0
        state['vi'] = - math.acos(
                        1 / eccentricity_in * (semimajoraxis_in * (1 - eccentricity_in ** 2) / r - 1))
        if args.montecarlo == True:
            from physical_models.MonteCarlo_perturbations import monte_carlo_true_anomaly
            state = monte_carlo_true_anomaly(state,args)
            apoapsis = state['Apoapsis']
            periapsis = state['Periapsis']

    # print('Apoapsis Radius: km', apoapsis / 10 ** 3, 'Periapsis Altitude: km', periapsis / 10 ** 3)


    ## Initial Model Definition
    # Body
    if args.body_shape == 'Spacecraft':
        def body ():
            Mass = mass  # kg
            length_SA = length_sp
            height_SA = height_sp
            Area_SA = length_SA*height_SA  # m^2
            length_SC = length_ody
            height_SC = height_ody
            Area_SC = length_ody*height_ody  # m^2
            Area_tot = Area_SA+Area_SC
            b = cnf.model.body(Mass, length_SA, height_SA, Area_SA, length_SC, height_SC, Area_SC, Area_tot)
            return b
    elif args.body_shape == 'Blunted Cone':
        def body():
            Mass = mass
            Delta = delta
            NoseRadius = nose_radius
            BaseRadius= base_radius
            Area_tot = math.pi*BaseRadius**2
            b = cnf.model.body(Mass, delta=Delta, NoseRadius=NoseRadius, BaseRadius=BaseRadius, Area_tot = Area_tot)
            return b
    b_class = body()


    # Initial Condition Entry
    def initialconditions ():
        a = semimajoraxis_in        #semi-major axis
        e = eccentricity_in         #eccentricity
        i = np.radians(state['Inclination'])                       #inclination # use 89.99 for 90
        OMEGA = np.radians(state['OMEGA'])                   #longitude of the ascending node
        omega = np.radians(state['omega'])                   #argument of the periapsis
        vi = state['vi']  # deg
        m0 = mass  # kg
        year = args.year
        month = args.month
        day = args.day
        hour = args.hours
        min = args.minutes
        second = args.secs
        time_rot = args.planettime
        ic = cnf.model.initialcondition(a , e , i , OMEGA , omega, vi, m0, year, month, day, hour, min, second, time_rot)
        return ic
    ic_class = initialconditions()

    # Aerodynamics
    def aerodynamics ():
        delta = math.radians(0)
        aoa = math.radians(args.angle_of_attack)
        thermal_accomodation_factor = args.thermal_accomodation_factor #0.2
        reflection_coefficient = args.reflection_coefficient # for diffuse reflection a_f =0, for specular reflection a_f = 1
        thermal_contact = 0 # for thermal perfect contact between rear and frontal area t_f = 0, for thermal insulation between rear and frontal area t_f = 1,
        heat_rate_limit = args.max_heat_rate #W/cm^2
        heat_load_limit = args.max_heat_load
        a = cnf.model.aerodynamics(delta, aoa, thermal_accomodation_factor, reflection_coefficient, thermal_contact, heat_rate_limit,heat_load_limit)
        return a
    a_class = aerodynamics()

    # Engine
    args.phi = math.radians(args.phi)
    def engine ():
        phi = args.phi  # deg # Thrust angle, angle between D and T, if T same direction of D = phi = 0
        g_e = 9.81  # m/s
        T = args.thrust  # N
        Isp = 200  # s
        e = cnf.model.engine(phi , g_e , T , Isp)
        return e
    e_class = engine()


    def model():
        body = b_class
        planet = p_class
        initialcondition = ic_class
        aerodynamics = a_class
        engine = e_class
        m = cnf.model(body, planet, initialcondition, aerodynamics, engine)
        return m
    m = model()

    #Initialization  # Reset all the config index for new simulation
    cnf.count_aerobraking = 0
    cnf.count_overcome_hr = 0
    cnf.save_index_heat = 0
    cnf.index_propellant_mass = 1
    cnf.counter_random = 0


    ##############################################
    # RUN SIMULATION                             #
    cnf.heat_rate_limit = args.max_heat_rate
    t = timecode.time()
    aerobraking(ip, m, args)
    elapsed = timecode.time() - t
    #############################################

    if args.print_res:
        print('rho: ',max(cnf.solution.physical_properties.rho), 'kg/m^3')
        print('heat rate: ',max(cnf.solution.performance.heat_rate), 'W/cm^2')


    # Save results
    if save_res == 1:
        now = datetime.now()
        if args.filename == 1:

            if args.montecarlo == True:
                folder_name = args.simulation_filename[0:args.simulation_filename.find('_nMC')]
            else:
                folder_name = args.simulation_filename
            name = args.directory_results + folder_name + '/' + args.simulation_filename
            filename = name + '.csv'

        else:
            name = args.directory_results + '/Sim' + str(args.MarsGram_version)
            # name = args.directory_results + str(date.today()) + '/' + str(now.strftime("%H_%M_%S"))

            filename = name + '.csv'

        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        save_cvs(filename,args)
    if args.print_res:
        print('Elapsed Time :', elapsed)

    if args.plot == True:
        from utils.plots import plots

        plots(state,m,name,args)









