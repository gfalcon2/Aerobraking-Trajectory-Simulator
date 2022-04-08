#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

import math


def ic_calculation_rptoae(planet, gamma, v, args):
    if args.drag_passage:
        h_0 = args.EI*1e3
    elif args.body_shape == 'Blunted Cone':
        h_0 = args.EI*1e3
    r = planet.Rp_e + h_0 # Drag passage always start and end at EI km of altitude
    a = planet.mu /((2.0*planet.mu/r)-v**2)
    h = r*v*math.cos(math.radians(gamma))
    p = (h**2)/planet.mu
    e = math.sqrt(1-p/a)
    ra = a*(1+e)
    hp = ((a*(1-e))- planet.Rp_e)  #(I'm considering a spherical planet)
    if hp < 0.0:
        print('WARNING AT initial_cond_calc: ALTITUDE PERIAPSIS < 0!')

    return ra, hp
