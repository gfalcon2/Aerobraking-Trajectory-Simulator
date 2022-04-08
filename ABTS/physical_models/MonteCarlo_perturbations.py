#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
## Note: HOW THE MONTECARLO PERTURBATIONS ARE SET UP
## the seed is the monte carlo number, the random number generator generates config.counter_random number, which increase of +1 every time the generator is called.
## the generates output is always the last number of the list. When config.counter_random number reaches 100, it is re-initialized to 1 to avoid storage problem and slowing down the function.
## Changing the config.counter_random number results in changing the random number generates, otherwise this would be the same everytime, biasing the results.
import numpy as np
from numpy import random
import config
def monte_carlo_aerodynamics(CL_body, CD_body,args): # Called multiple times during the simulation - for each time step
    uncertainty_CD,uncertainty_CL = args.CD_dispersion/100, args.CL_dispersion/100
    CD_body += uniform_distribution(CD_body*uncertainty_CD,number = 1)
    CL_body += uniform_distribution(CL_body*uncertainty_CL,number = 2)
    return CL_body, CD_body

def monte_carlo_density(density, args): # Called multiple times during the simulation - for each time step
    # uncertainty_rho = args.CD_dispersion/100, args.CL_dispersion/100
    random.seed(int(config.index_MonteCarlo))
    # rho = density
    # print()
    density = random.uniform(density * 0.5, density*2, 3)[-1]
    # print(rho,density,random.uniform(density * 0.7, density*1.5, 3)[-1])
    # pert = random.uniform(density,1.1)[-1]
    # rho = density*pert

    # print(pert)
    # CL_body += uniform_distribution(CL_body*uncertainty_CL,number = 2)
    return density

def monte_carlo_initial_conditions(state,args): # Called one time
    state['Apoapsis'] += uniform_distribution(args.ra_dispersion,number = 3)#r[0]
    state['Periapsis'] += uniform_distribution(args.rp_dispersion,number = 4)#r[1]
    state['Inclination'] += uniform_distribution(args.i_dispersion,number = 5)#s[0]
    state['OMEGA'] += uniform_distribution(args.raan_dispersion,number = 6)# s[1]
    state['omega'] += uniform_distribution(args.aop_dispersion,number = 7)#s[2]
    state['vi'] += uniform_distribution(args.vi_dispersion,number = 8)
    return state

def monte_carlo_true_anomaly(state,args): # Called one time
    state['vi'] += uniform_distribution(args.vi_dispersion,number = 9)#s[2]
    return state

def monte_carlo_guidance_closedform(state,args):  # Called one time
    state['ra'] += uniform_distribution(args.ra_dispersion_gnc,number = 10 )#r[0]
    state['rp'] += uniform_distribution(args.rp_dispersion_gnc,number = 11)#r[1]
    state['i'] += uniform_distribution(np.radians(args.i_dispersion_gnc),number = 12)#s[0]
    state['OMEGA'] += uniform_distribution(np.radians(args.raan_dispersion_gnc),number = 13)# s[1]
    state['omega'] += uniform_distribution(np.radians(args.aop_dispersion_gnc),number = 14)#s[2]
    state['vi'] += uniform_distribution(np.radians(args.vi_dispersion_gnc),number = 15)
    return state

def monte_carlo_guidance_environment(rho,T,S,args):  # Called multiple times during the simulation - with the guidance frequency

    mean_rho,mean_T,mean_S = args.rho_mudispersion_gnc/100, args.T_mudispersion_gnc/100,args.S_mudispersion_gnc/100
    std_rho,std_T,std_S = args.rho_sigmadispersion_gnc/100, args.T_sigmadispersion_gnc/100,args.S_sigmadispersion_gnc/100

    rho += gaussian_distribution(rho*mean_rho,rho*std_rho,number=16)#)config.counter_random)
    T += gaussian_distribution(T*mean_T,T*std_T,number=17)#config.counter_random)
    S += gaussian_distribution(S*mean_S,S*std_S,number=18)#config.counter_random)
    return rho, T, S

def uniform_distribution(dispersion, number = None):
    random.seed(int(config.index_MonteCarlo))
    return random.uniform(-dispersion,dispersion,number)[-1]

def gaussian_distribution(mu, sigma, number = None):
    config.counter_random += 1
    random.seed(int(config.counter_random))
    number = (random.randint(1,500))
    random.seed(int(config.index_MonteCarlo))
    return random.normal(mu , sigma , number)[-1]


