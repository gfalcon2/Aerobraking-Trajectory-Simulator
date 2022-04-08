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
from utils.Odyssey_maneuver_plan import Odyssey_firing_plan
from utils.save_results import clean_results
from physical_models.Propulsive_maneuvers import propulsion_ic_calcs

def aerobraking(ip, m, args):
    initial_state = m.initialcondition
    FinalState = True
    continue_campaign = True
    numberofpassage = 0    # Initialization
    cnf.count_numberofpassage = 0
    clean_results() # Initialize Solution Results
    cnf.time_OP = 0
    cnf.time_IP = 0
    # Aerobraking Campaign
    while continue_campaign and (FinalState):
        cnf.index_Mars_Gram_call = 0
        cnf.firing_orbit = 0
        numberofpassage = 1 + numberofpassage
        if args.print_res:
            print("--> START PASSAGE #{}".format(numberofpassage))
        t = time.time()

        if args.Odyssey_sim:
            args = Odyssey_firing_plan(numberofpassage,args)
        if ip.tc == 1:
            # calculate how long before the apoapsis the maneuver should start
            if args.delta_v != 0.0:
                if ip.tc == 1:
                    initial_state = propulsion_ic_calcs(m, args, initial_state)
        elif ip.tc == 2:
            if round(math.degrees(args.phi)) == 180:
                print("DECELERATE DRAG FIRING!!")
            elif round(math.degrees(args.phi)) == 0:
                print("ACCELERATE DRAG FIRING!")

        if numberofpassage != 1:
            # Orbital Elements Result
            initial_state.a = cnf.solution.orientation.oe[0][-1]
            initial_state.e = cnf.solution.orientation.oe[1][-1]
            initial_state.i = cnf.solution.orientation.oe[2][-1]
            initial_state.OMEGA = cnf.solution.orientation.oe[3][-1]
            initial_state.omega = cnf.solution.orientation.oe[4][-1]
            initial_state.m = cnf.solution.performance.mass[-1]
            initial_state.vi = cnf.solution.orientation.oe[5][-1]

            m.initialcondition.year = int(cnf.solution.orientation.year[-1])
            m.initialcondition.month = int(cnf.solution.orientation.month[-1])
            m.initialcondition.day = int(cnf.solution.orientation.day[-1])
            m.initialcondition.hour = int(cnf.solution.orientation.hour[-1])
            m.initialcondition.min = int(cnf.solution.orientation.min[-1])
            m.initialcondition.second = int(cnf.solution.orientation.second[-1])

            if (args.drag_passage or args.body_shape == 'Blunted Cone') and continue_campaign:
                r = m.planet.Rp_e + args.EI*10**3
                initial_state.vi = - math.acos(1 / initial_state.e * (initial_state.a * (1 - initial_state.e ** 2) / r - 1))

        # Run Simulation
        continue_campaign = asim(ip, m, initial_state, numberofpassage, args)

        r_a = cnf.solution.orientation.oe[0][-1] * (1 + cnf.solution.orientation.oe[1][-1])
        r_p = cnf.solution.orientation.oe[0][-1]* (1 - cnf.solution.orientation.oe[1][-1])
        elapsed = time.time() - t

        if args.print_res:
            print('Computational Time {}s'.format(elapsed))
            print("--> PASSAGE #{} COMPLETE".format(numberofpassage))

        if args.number_of_orbits == numberofpassage:
            continue_campaign = False

        if r_a <= args.final_apoapsis:
            FinalState = False
            print('REACHED FINALSTATE! R_a = ',r_a*1e-3,'km')
            print('Thermal Limit overcomed totally', cnf.count_overcome_hr, 'times!')

        if r_p - m.planet.Rp_e>= 180*1e3:
            FinalState = False
            print('PERIAPSIS TOO HIGH, FINAL STATE UNREACHABLE! R_a = ', r_p * 1e-3, 'km')

        print(' ')

    closed_form(args, m)

