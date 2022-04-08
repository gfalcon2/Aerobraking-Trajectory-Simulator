#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import config
import math
# F = N = kg*m/s^2 = T*t/mass (consider mass constant in the process)
# v = m/s = kg*m/s^2 * s/kg = T *t /m
# t = v*m/T
def Odyssey_firing_plan(numberofpassage,args):
    if numberofpassage == 7:
        args.delta_v = 0.14
        args.phi = math.radians(180)
    elif numberofpassage == 14: #13
        args.delta_v =  0.15
        args.phi = 0
    elif numberofpassage == 26:
        args.delta_v =  0.1
        args.phi = 0
    elif numberofpassage == 30:
        args.delta_v =  0.1
        args.phi = 0
    elif numberofpassage == 35:
        args.delta_v = 0.2
        args.phi = 0
    elif numberofpassage == 47:
        args.delta_v = 0.2
        args.phi = 0
    elif numberofpassage == 55:
        args.delta_v =  0.3
        args.phi = math.radians(180)
    elif numberofpassage == 69:
        args.delta_v = 0.15
        args.phi = 0
    elif numberofpassage == 72:
        args.delta_v = 0.15
        args.phi = 0
    elif numberofpassage == 80:
        args.delta_v = 0.15
        args.phi = 0
    elif numberofpassage == 87:
        args.delta_v = 0.14/2
        args.phi = math.radians(180)
    elif numberofpassage == 110:
        args.delta_v = 0.14/2
        args.phi = math.radians(180)
    elif numberofpassage == 127:
        args.delta_v = 1.0/3
        args.phi = math.radians(180)
    elif numberofpassage == 160:
        args.delta_v = 0.84/3
        args.phi = math.radians(180)
    elif numberofpassage == 178:
        args.delta_v = 0.6/3
        args.phi = math.radians(180)
    elif numberofpassage == 194:
        args.delta_v = 0.84/3
        args.phi = math.radians(180)
    elif numberofpassage == 210:
        args.delta_v = 0.6/3
        args.phi = math.radians(180)
    elif numberofpassage == 222:
        args.delta_v = 0.6/3
        args.phi = math.radians(180)
    elif numberofpassage == 238:
        args.delta_v = 1.2/3
        args.phi = math.radians(180)
    elif numberofpassage == 251:
        args.delta_v = 1.0/3
        args.phi = math.radians(180)
    elif numberofpassage == 263:
        args.delta_v = 1./3
        args.phi = math.radians(180)
    elif numberofpassage == 274:
        args.delta_v = 1.2/3
        args.phi = math.radians(180)
    elif numberofpassage == 287:
        args.delta_v = 1/3
        args.phi = math.radians(180)
    elif numberofpassage == 299:
        args.delta_v = 1/3
        args.phi = math.radians(180)
    elif numberofpassage == 311:
        args.delta_v = 1.2/3
        args.phi = math.radians(180)
    else:
        args.delta_v = 0.0
        args.phi = 0.0

    # if args.phi == math.radians(0) and args.delta_v != 0.0:
    #     print("LOWER MANEUVER!")
    # elif args.phi == math.radians(180) and args.delta_v != 0.0:
    #     print("RAISE MANEUVER!")

    return args
#orbit_6 = raise 0.18m/s
#orbit 13 = lower 0.18 m/s
#orbit 22 lower 0.1
#orbit 31 lower 0.1
#orbit 36 lower 0.2
# orbit 48 lower 0.2
# orbit 56 raise 0.3
#  orbit 69 lower 0.18
# orbit 71 lower  0.18
# orbit 81 lower 0.18
# orbit 86 raise 0.18
#orbit 111 raise  0.18
#orbit 126 raise 1
#orbit 161 raise 0.85
#orbit 176 raise 0.6
#orbit 183 raise 0.85
#orbit 211 raise 0.6
# orbit 221 raise 0.6