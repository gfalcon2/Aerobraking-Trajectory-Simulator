#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import config
import numpy as np
import csv
import os

def MonteCarlo_setting(args):
    MC = {}
    MC['N passages'] , MC['Duration'] , MC['Median Heat'] , MC['Periapsis min'] , MC[
        'Periapsis max'] , count = [] , [] , [] , [] , [] , 0
    if args.montecarlo == False:
        args.montecarlo_size = 2
        args.initial_montecarlo_number = 1
    return MC, count, args

def MonteCarlo_setting_passage(mc_index,args):
    config.counter_random = 0
    heat_passage = []
    if args.print_res:
        print('--> MC number' , mc_index)
    # if commented: filename = current time (comment if you don't need)
    args.simulation_filename = args.simulation_filename + '_nMC={}'.format(mc_index)


    # Initialization
    config.altitudeperiapsis , config.max_heatrate = [] , []
    config.counter_random = mc_index
    config.index_MonteCarlo = mc_index
    return args

def MonteCarlo_append(MC, args , count):
    # append results Montecarlo
    config.index_MonteCarlo += 1

    # Save results
    MC['N passages'].append(config.solution.orientation.numberofpassage[-1])
    MC['Duration'].append(config.solution.orientation.time[-1])
    MC['Median Heat'].append(np.median(config.max_heatrate))
    MC['Periapsis min'].append(min(config.altitudeperiapsis))
    MC['Periapsis max'].append(max(config.altitudeperiapsis))
    heat_rate_max = max(config.solution.performance.heat_rate)
    if heat_rate_max > args.max_heat_rate:
        count += 1
    if args.print_res:
        print('--> Count =' , count)


def MonteCarlo_save(args,state,MC):
    if args.montecarlo_analysis:
        # save results montecarlo, comparison, only if analysis is in process
        folder_name = args.simulation_filename[0:args.simulation_filename.find('_nMC')]

        name = args.directory_results+ folder_name + '/MC_results_control={}'.format(args.control_mode)+'_ra={}'.format(int(state['Apoapsis']/1e3))+'_rp={0:.1f}'.format(state['Periapsis'])+'_hl={}'.format(args.max_heat_rate)+'.csv'
        filename = name + '.csv'

        os.makedirs(os.path.dirname(filename) , exist_ok=True)
        with open(filename , "w") as f:
            writer = csv.writer(f , delimiter=',' , quotechar='"' , quoting=csv.QUOTE_MINIMAL)
            writer.writerow(range(args.montecarlo_size))
            writer.writerow(MC['N passages'])
            writer.writerow(MC['Duration'])
            writer.writerow(MC['Median Heat'])
            writer.writerow(MC['Periapsis min'])
            writer.writerow(MC['Periapsis max'])

