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
import config

cpdef run_orbitalelements(object args):
    #initialization
    # cdef int t = timecode.time()
    cdef list apoapsis, periapsis_alt
    cdef double inclination, OMEGA, omega, periapsis_item, apoapsis_item
    cdef dict state, MC
    cdef int count, mc_index

    apoapsis , periapsis_alt , inclination ,OMEGA, omega = list(
    range(int(args.ra_initial_a) , int(args.ra_initial_b) , int(args.ra_step))) , \
    list(range(int(args.hp_initial_a) , int(args.hp_initial_b) ,
    int(args.hp_step))) , args.inclination , args.OMEGA, args.omega  # 93.6 , 0 28523.95142557378  # with values of omega not equal to 0, recalculate periapsis altitude due to the oblateness of the planet
    cdef double final_apoapsis = args.final_apoapsis

    for periapsis_item in periapsis_alt:
        for apoapsis_item in apoapsis:
            if args.print_res:
                print('Apoapsis Radius: km' , apoapsis_item / 10 ** 3 , 'Periapsis Altitude: km' , periapsis_item / 10 ** 3)
            state = {}
            MC,count,args = MonteCarlo_setting(args)
            for mc_index in range(args.initial_montecarlo_number , args.montecarlo_size):
                state['Apoapsis'] , state['Periapsis'] , state['Inclination'] , state['OMEGA'] , state[
                    'omega'] , state[
                    'Final SMA'] = apoapsis_item , float(
                    periapsis_item)*1e-3 , inclination , OMEGA , omega , final_apoapsis

                args.simulation_filename = 'Results_ctrl={}'.format(args.control_mode) + '_ra={}'.format(
                    int(apoapsis_item*1e-3)) + '_rp={0:.1f}'.format(
                    float(periapsis_item*1e-3)) + '_hl={0:.3f}'.format(
                    args.max_heat_rate) + '_{0:.1f}'.format(args.angle_of_attack)+'deg'

                if args.montecarlo == True:
                    args = MonteCarlo_setting_passage(mc_index , args)

                aeroraking_campaign(args,state)
                MonteCarlo_append(MC , args , count)

            if args.montecarlo == True:
                MonteCarlo_save(args, state,MC)

    # print('--> COMPUTATIONAL TIME = ' , timecode.time() - t)


cpdef run_vgamma(object args):
    #initialization
    # cdef int t = timecode.time()
    cdef list gamma_0, v_0
    cdef double inclination, OMEGA, omega, gamma, v, apoapsis, periapsis_alt
    cdef dict state, MC
    cdef int count, mc_index
    cdef object planet

    gamma_0 , v_0 , inclination , OMEGA, omega = list(range(int(args.gamma_initial_a*100) , int(args.gamma_initial_b*100) , int(args.gamma_step*100))) , \
    list(range(int(args.v_initial_a) , int(args.v_initial_b) ,
    int(args.v_step))) , args.inclination ,args.OMEGA, args.omega  # 93.6 , 0 28523.95142557378  # with values of omega not equal to 0, recalculate periapsis altitude due to the oblateness of the planet
    cdef double final_apoapsis = args.final_apoapsis

    for gamma in gamma_0:
        gamma = -gamma / 100
        for v in v_0:
            state={}
            planet = planet_data(args.planet)
            apoapsis , periapsis_alt = ic_calculation_rptoae(planet , gamma , v, args)
            if args.print_res:
                print('Velocity:', v,' m/s, Flight-Path Angle: ',gamma,'deg')

            MC,count,args = MonteCarlo_setting(args)

            for mc_index in range(args.initial_montecarlo_number , args.montecarlo_size):
                state['Apoapsis'] , state['Periapsis'] , state['Inclination'] , state['OMEGA'] , state[
                    'omega'] , state[
                    'Final SMA'] = apoapsis , float(
                    periapsis_alt) *1e-3 , inclination , OMEGA , omega , final_apoapsis

                args.simulation_filename = 'Results_ctrl={}'.format(args.control_mode) + '_v={0:.2f}'.format(
                    int(v)) + '_gamma={0:.2f}'.format(abs(gamma)) + '_{0:.1f}'.format(args.angle_of_attack)+'deg'
                if args.montecarlo == True:
                    args = MonteCarlo_setting_passage(mc_index , args)

                aeroraking_campaign(args,state)
                MonteCarlo_append(MC , args , count)

            if args.montecarlo == True:
                MonteCarlo_save(args, state,MC)


    # print('--> COMPUTATIONAL TIME = ' , timecode.time() - t)

cpdef run_analysis(object args):

    args = def_miss(args)
    if args.initial_condition_type==1 and (args.drag_passage or args.body_shape == 'Blunted Cone'):
        run_vgamma(args)
    else:
        run_orbitalelements(args)

    if args.passresults:
        return config.solution
    else:
        return bool(1)
