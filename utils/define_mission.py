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
        args.number_of_orbit = 1000

    if args.body_shape == 'Spacecraft':
        args.aerodynamic_model = 'Mach-dependent'
        args.thermal_model = 'Maxwellian Heat Transfer'
    elif args.body_shape == 'Blunted Cone':
        args.aerodynamic_model = 'No-Ballistic flight with axial coefficient'
        args.thermal_model = 'Convective'

    if args.thrust_control == 'None':
        args.thrust = 0
    elif args.thrust_control == 'Aerobraking Maneuver':
        args.thrust = args.thrust
    elif args.thrust_control == 'Drag Passage Firing':
        args.thrust = args.thrust




    return args