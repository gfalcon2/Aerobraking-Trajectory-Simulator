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
        delta = 70 #deg
        nose_radius = 0.73 #m
        base_radius = 0.75 #m


    apoapsis = state['Apoapsis']# 5000*10**3#7523.95142557378*10**3#28523.95142557378*10**3 # From phase 1 30371.25610549734*10**3#43090.01227466145 *10**3 #30371.25610549734*10**3#25000 *10**3 # from Phase 2 28523.95142557378*10**3 km

    state['Periapsis'] = p_class.Rp_e + (state['Periapsis'])*10**3 #98
    if args.montecarlo == True:
        from physical_models.MonteCarlo_perturbations import monte_carlo_initial_conditions
        state = monte_carlo_initial_conditions(state,args)



    # Initial Condition Calcs
    semimajoraxis_in = (state['Apoapsis'] + state['Periapsis'])/2
    eccentricity_in = (state['Apoapsis'] - state['Periapsis'])/ (state['Apoapsis'] + state['Periapsis'])
    state['vi'] = np.radians(180.0001)

    # Initial Condition
    if args.drag_passage == True:
        #if (simulation['IC v-r'] == True):  # The initial conditions in orbital elements need to be calculated:
        r = p_class.Rp_e + 160 * 10 ** 3
        state['vi'] = - math.acos(
                        1 / eccentricity_in * (semimajoraxis_in * (1 - eccentricity_in ** 2) / r - 1))

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
            b = cnf.model.body(Mass, length_SA, height_SA, Area_SA, length_SC, height_SC, Area_SC)
            return b
    elif args.body_shape == 'Blunted Cone':
        def body():
            Mass = mass
            Delta = delta
            NoseRadius = nose_radius
            BaseRadius= base_radius
            b = cnf.model.body(Mass, delta=Delta, NoseRadius=NoseRadius, BaseRadius=BaseRadius)
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
        ic = cnf.model.initialcondition(a , e , i , OMEGA , omega, vi, m0, year, month, day, hour, min, second)
        return ic
    ic_class = initialconditions()

    # Aerodynamics
    def aerodynamics ():
        delta = math.radians(0)
        aoa = math.radians(args.angle_of_attack)
        thermal_accomodation_factor = args.thermal_accomodation_factor #0.2
        accomodation_factor = args.accomodation_factor # for diffuse reflection a_f =0, for specular reflection a_f = 1
        thermal_contact = 0 # for thermal perfect contact between rear and frontal area t_f = 0, for thermal insulation between rear and frontal area t_f = 1,
        thermal_limit = args.max_heat_rate #W/cm^2
        a = cnf.model.aerodynamics(delta, aoa, thermal_accomodation_factor, accomodation_factor, thermal_contact, thermal_limit)
        return a
    a_class = aerodynamics()

    # Engine
    def engine ():
        phi = math.radians(args.phi)  # deg # Thrust angle, angle between D and T, if T same direction of D = phi = 0
        g_e = 9.81  # m/s
        T = args.thrust  # N
        Isp = 300  # s
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
                folder_name= args.simulation_filename
            name = args.directory_results + folder_name + '/' + args.simulation_filename
            filename = name + '.csv'

        else:
            name = args.directory_results + str(date.today()) + '/' + str(now.strftime("%H_%M_%S"))

            filename = name + '.csv'

        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        save_cvs(filename,args)
    if args.print_res:
        print('Elapsed Time :', elapsed)

    if args.plot == True:
        from utils.plots import plots

        plots(state,m,name,args)








