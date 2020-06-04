#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

import config as cnf
import math
from simulation.Complete_passage import asim
from utils.Ref_system_conf import OE
import time
import numpy as np
from utils.closed_form_solution import closed_form

def aerobraking(ip, m, args):
    initial_state = m.initialcondition
    FinalState = True
    continue_campaign = True
    numberofpassage = 0
    time_0 = 0

    # Aerobraking Campaign
    #while numberofpassage<1:
    while continue_campaign and (FinalState):
        cnf.index_Mars_Gram_call = 0
        cnf.firing_orbit = 0
        numberofpassage = 1 + numberofpassage
        print("--> START PASSAGE #{}".format(numberofpassage))
        t = time.time()

        # Run Simulation
        continue_campaign = asim(ip, m, time_0, initial_state, numberofpassage, args)

        # Orbital Elements Result
        a = cnf.solution.orientation.oe[0][-1]
        e = cnf.solution.orientation.oe[1][-1]
        i = cnf.solution.orientation.oe[2][-1]
        OMEGA = cnf.solution.orientation.oe[3][-1]
        omega = cnf.solution.orientation.oe[4][-1]
        mass = cnf.solution.performance.mass[-1]
        vi = math.pi+0.01

        # New Initial State
        initial_state = OE(a, e, i, OMEGA, omega, vi, mass)

        v = [cnf.solution.orientation.vel_ii[0][-1], cnf.solution.orientation.vel_ii[1][-1], cnf.solution.orientation.vel_ii[2][-1]]
        r = [cnf.solution.orientation.pos_ii[0][-1], cnf.solution.orientation.pos_ii[1][-1], cnf.solution.orientation.pos_ii[2][-1]]
        time_0 = cnf.solution.orientation.time[-1]
        r_a = a*(1+e)
        elapsed = time.time() - t


        #print('Mars Gram called {} times'.format(cnf.index_Mars_Gram_call))
        print('Computational Time {}s'.format(elapsed))
        print("--> PASSAGE #{} COMPLETE".format(numberofpassage))
        #print('Energy', ((np.inner(v, v)) / 2 - m.planet.mu / (np.inner(r, r)) ** 0.5))

        if args.number_of_orbits == numberofpassage:
            continue_campaign = False

        if r_a <= args.final_apoapsis:
            FinalState = False
            print('Reached Final State! R_a = ',r_a,'km')
            print('Thermal Limit overcomed totally', cnf.count_overcome_hr, 'times!')
        print(' ')

    closed_form(args, m)
