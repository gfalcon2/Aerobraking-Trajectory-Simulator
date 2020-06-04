#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

from simulation.Set_and_run import *
import time as timecode
from utils.initial_cond_calc import *
from utils.MonteCarlo_set import *
import os
import sys
import io
from utils.define_mission import def_miss

def run_orbitalelements(args):
    t = timecode.time()
    apoapsis , periapsis_alt , inclination , omega = list(
    range(int(args.ra_initial_a) , int(args.ra_initial_b) , int(args.ra_step))) , \
    list(range(int(args.hp_initial_a) , int(args.hp_initial_b) ,
    int(args.hp_step))) , args.inclination , args.omega  # 93.6 , 0 28523.95142557378  # with values of omega not equal to 0, recalculate periapsis altitude due to the oblateness of the planet
    final_apoapsis = args.final_apoapsis


    for periapsis_item in periapsis_alt:
        for apoapsis_item in apoapsis:
            print('Apoapsis Radius: km' , apoapsis_item / 10 ** 3 , 'Periapsis Altitude: km' , periapsis_item / 10 ** 3)
            state = {}

            if args.montecarlo == False:
                args.montecarlo_size = 1

            MC,count = MonteCarlo_setting()

            for mc_index in range(0 , args.montecarlo_size):
                state['Apoapsis'] , state['Periapsis'] , state['Inclination'] , state['OMEGA'] , state[
                    'omega'] , state[
                    'Final SMA'] = apoapsis_item , float(
                    periapsis_item) / 10 ** 3 , inclination , omega , omega , final_apoapsis

                args.simulation_filename = 'Results_ctrl={}'.format(args.control_mode) + '_ra={}'.format(
                    int(apoapsis_item / 10 ** 3)) + '_rp={0:.1f}'.format(
                    float(periapsis_item / 10 ** 3)) + '_hl={}'.format(
                    args.max_heat_rate) + '_0deg'

                if args.montecarlo == True:
                    args = MonteCarlo_setting_passage(mc_index , args)

                aeroraking_campaign(args,state)
                MonteCarlo_append(MC , args , count)

            if args.montecarlo == True:
                MonteCarlo_save(args, state,MC)

    print('--> COMPUTATIONAL TIME = ' , timecode.time() - t)


def run_vgamma(args):
    t = timecode.time()
    gamma_0 , v_0 , inclination , omega = list(range(int(args.gamma_initial_a*100) , int(args.gamma_initial_b*100) , int(args.gamma_step*100))) , \
    list(range((args.v_initial_a) , (args.v_initial_b) ,
    (args.v_step))) , args.inclination , args.omega  # 93.6 , 0 28523.95142557378  # with values of omega not equal to 0, recalculate periapsis altitude due to the oblateness of the planet
    final_apoapsis = args.final_apoapsis

    for gamma in gamma_0:
        gamma = -gamma / 100
        for v in v_0:
            state={}
            planet = planet_data(args.planet)
            apoapsis , periapsis_alt = ic_calculation_rptoae(planet , gamma , v)
            inclination , omega = args.inclination, args.omega

            print('Apoapsis Radius: km' , apoapsis / 10 ** 3 , 'Periapsis Altitude: km' , periapsis_alt / 10 ** 3)
            print('V', v,'Gamma',gamma)


            if args.montecarlo == False:
                args.montecarlo_size = 1

            MC,count = MonteCarlo_setting()

            for mc_index in range(0 , args.montecarlo_size):

                # Orbital Element
                state['Apoapsis'] , state['Periapsis'] , state['Inclination'] , state['OMEGA'] , state[
                    'omega'] , state[
                    'Final SMA'] = apoapsis , float(
                    periapsis_alt) / 10 ** 3 , inclination , omega , omega , final_apoapsis

                args.simulation_filename = 'Results_ctrl={}'.format(args.control_mode) + '_v={}'.format(
                    int(v)) + '_gamma={}'.format(abs(gamma)) + '_90deg'
                args = MonteCarlo_setting_passage(mc_index , args)

                aeroraking_campaign(args,state)
                MonteCarlo_append(MC , args , count)

            if args.montecarlo == True:
                MonteCarlo_save(args, state,MC)


    print('--> COMPUTATIONAL TIME = ' , timecode.time() - t)

def run_analysis(args):
    args = def_miss(args)
    if (args.drag_passage == True) and (args.initial_condition_type== 1):
        run_vgamma(args)
    else:
        run_orbitalelements(args)