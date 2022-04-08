#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

def def_miss(args):

    if args.type_of_mission == 'Drag Passage':
        args.drag_passage = 1
        args.number_of_orbits = 1
    elif args.type_of_mission == 'Orbits':
        args.drag_passage = 0
        args.number_of_orbits = args.number_of_orbits
    elif args.type_of_mission == 'Aerobraking Campaign':
        args.drag_passage = 0
        args.number_of_orbits = 1000

    if args.body_shape == 'Spacecraft':
        if args.aerodynamic_model == 'No-Ballistic flight with axial coefficient':
            args.aerodynamic_model = 'Mach-dependent'
            print('--AERODYNAMIC MODEL CHANGED TO: Mach-dependent - Specific for a flat-plate--')
        if args.thermal_model != 'Maxwellian Heat Transfer':
            args.thermal_model = 'Maxwellian Heat Transfer'
            print('--THERMAL MODEL CHANGED TO: Maxwellian Heat Transfer - Specific for a flat-plate--')
    elif args.body_shape == 'Blunted Cone':
        if args.aerodynamic_model == 'Mach-dependent':
            print('--AERODYNAMIC MODEL CHANGED TO: No-Ballistic flight with axial coefficient - Specific for a Blunted-Cone--')
            args.aerodynamic_model = 'No-Ballistic flight with axial coefficient'
        if args.thermal_model != 'Convective and Radiative':
            args.thermal_model = 'Convective and Radiative'
            print('--THERMAL MODEL CHANGED TO: Convective and Radiative - Specific for a Blunted-Cone--')
        if args.control_mode != 0:
            print('--ARTICULATED SOLAR PANELS GUIDANCE NOT ALLOWED FOR BLUNTED CONE--')
            args.control_mode = 0

    if args.thrust_control == 'None':
        args.thrust = args.thrust
        args.delta_v = 0
    elif args.thrust_control == 'Aerobraking Maneuver':
        args.thrust = args.thrust
        args.delta_v = args.delta_v
        if args.type_of_mission == 'Drag Passage':
            args.thrust_control = 'None'
    elif args.thrust_control == 'Drag Passage Firing':
        args.thrust = args.thrust
        args.delta_v = args.delta_v
    if args.thrust_control != 'None' and args.thrust == 0:
        print('THRUST MODIFIED TO 0.1 N')
        args.thrust = 0.1

    if args.machine == 'Aero':
        if not args.directory_results:
            args.directory_results = '/Users/giusyfalcone/Aerobraking_SA_project_results/'
        args.directory_MarsGram = '/Users/giusyfalcone/MarsGram/'
    elif args.machine == 'Cluster':
        if not args.directory_results:
            args.directory_results = '/home/gfalcon2/scratch/Aerobraking_SA_project_results/'
        args.directory_MarsGram = '/home/gfalcon2/scratch/Aerobraking/ABTS/MarsGram/'
    elif args.machine == 'Desktop_Home':
        if not args.directory_results:
            args.directory_results = '/Users/Giusy/Aerobraking_SA_project_results/'
        args.directory_MarsGram = '/Users/Giusy/MarsGram/'
    elif args.machine == 'Laptop':
        if not args.directory_results:
            args.directory_results = '/Users/Josephine/Aerobraking/Aerobraking_SA_project_results/'
        args.directory_MarsGram = '/Users/Josephine/Aerobraking/MarsGram/'

    if args.density_model == 'MarsGram' and args.planet != 1:
        args.density_model = 'Exponential'
        print('Density Model Changed with Exponential Law - GRAM only with Mars')

    if args.Odyssey_sim:
        args.control_mode = 0
        args.type_of_mission = 'Aerobraking Campaign'
        args.number_of_orbits = 1000
        args.planet = 'Mars'
        args.body_shape = 'Spacecraft'
        args.dry_mass = 411.0
        args.prop_mass = 50.0
        args.angle_of_attack = 90.0
        args.initial_condition_type = 0
        args.thrust_control = 'Aerobraking Maneuver'
        args.ra_initial_a = 28523.95 * 10 ** 3
        args.ra_initial_b = 30000.00 * 10 ** 3
        args.ra_step = 1e12
        if args.gravity_model == 'Inverse Squared':
            args.hp_initial_a = 108600
        else:
            args.hp_initial_a = 86000#98600 #84280#
        args.hp_initial_b = 120000
        args.hp_step = 1e12
        args.inclination = 93.6
        args.omega = 109.0
        args.OMEGA = 114.0
        args.year = 2001
        args.month = 11
        args.day = 6
        args.final_apoapsis = 3900 * 10 ** 3
        args.montecarlo = 0
        args.drag_passage = 0
    return args
