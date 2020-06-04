#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

from numpy import random
import config
def monte_carlo_aerodynamics(CL_body, CD_body,args):
    random.seed(int(config.index_MonteCarlo))
    uncertainty_CD,uncertainty_CL = args.CD_dispersion/100, args.CL_dispersion/100
    CD_body += uniform_distribution(CD_body*uncertainty_CD)
    CL_body += uniform_distribution(CL_body*uncertainty_CL)
    config.counter_random = 0
    return CL_body, CD_body

def monte_carlo_initial_conditions(state,args):
    state['Apoapsis'] += uniform_distribution(args.ra_dispersion)#r[0]
    state['Periapsis'] += uniform_distribution(args.rp_dispersion)#r[1]
    state['Inclination'] += uniform_distribution(args.i_dispersion)#s[0]
    state['OMEGA'] += uniform_distribution(args.raan_dispersion)# s[1]
    state['omega'] += uniform_distribution(args.aop_dispersion)#s[2]
    return state

def uniform_distribution(dispersion):
    random.seed(int(config.index_MonteCarlo))
    if config.counter_random >= 10:
        config.counter_random = 0
    config.counter_random += 1
    return random.uniform(-dispersion,dispersion,config.counter_random)[-1]