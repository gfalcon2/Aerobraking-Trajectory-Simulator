#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

import config
import numpy as np
# from datetime import datetime

def save_results(time,ratio):
    # startTime = datetime.now()
    initial_time = 0
    if len(config.solution.orientation.time) == 0:
        config.prev_step_integrator = 0.0
        initial_time = config.prev_step_integrator
        config.initial_time_save = 0.0
    else:
        if config.counter_integrator == 1:
            config.prev_step_integrator = 0.0
            initial_time = config.prev_step_integrator
            config.initial_time_save = config.solution.orientation.time[-1]

    n_variable_to_save = len(config.solution_intermediate[0])-1 # number defined in the Complete_passage.py
    range_time = [item[0] for item in config.solution_intermediate]
    results = np.empty([n_variable_to_save ,1])
    t = []
    i = 0
    index_prev = 0
    for true_time in time: # total solution time
        if i % ratio == 0:
            try:
                index = range_time[index_prev:-1].index(true_time)
                t.append(range_time[index+index_prev]+initial_time)
                range_solution = np.reshape(config.solution_intermediate[index+index_prev] , (n_variable_to_save + 1 , 1))
                if len(t) == 1:
                    results = range_solution[1:]
                else:
                    results = np.append(results , range_solution[1:] , axis=1)
                index_prev = index # look for only after the index already used
            except:
                continue
                # if true_time == time[-1]:
                #     from utils.misc import find_neighbour
                #     index_a, index_b, dist_a, dist_b = find_neighbour(true_time,range_time[index_prev:-1])
                #     t.append((range_time[index_a]+(range_time[index_b]-range_time[index_a])/(dist_b+abs(dist_a)) * abs(dist_a))+initial_time)
                #     sol_a = np.reshape(config.solution_intermediate[index_a] , (n_variable_to_save + 1 , 1))
                #     sol_b = np.reshape(config.solution_intermediate[index_b] , (n_variable_to_save + 1 , 1))
                #     range_solution = [sol_a[i]+(sol_b[i]-sol_a[i])/(dist_b+abs(dist_a)) * abs(dist_a) for i in range(len(sol_a))]
                #     results = np.append(results , range_solution[1:] , axis=1)
                # else:
                #     continue
        i += 1

    time_0 = t[-1]
    config.prev_step_integrator = time_0
    config.solution_intermediate = []

    t = [i+config.initial_time_save for i in t]


    ## SAVE RESULTS GLOBAL
    # Orientation
    config.solution.orientation.time.extend(t)
    config.solution.orientation.year.extend(results[0])
    config.solution.orientation.month.extend(results[1])
    config.solution.orientation.day.extend(results[2])
    config.solution.orientation.hour.extend(results[3])
    config.solution.orientation.min.extend(results[4])
    config.solution.orientation.second.extend(results[5])
    config.solution.orientation.numberofpassage.extend(results[6])
    config.solution.orientation.pos_ii[0].extend(results[7])
    config.solution.orientation.pos_ii[1].extend(results[8])
    config.solution.orientation.pos_ii[2].extend(results[9])
    config.solution.orientation.vel_ii[0].extend(results[10])
    config.solution.orientation.vel_ii[1].extend(results[11])
    config.solution.orientation.vel_ii[2].extend(results[12])
    config.solution.orientation.pos_ii_mag.extend(results[13])
    config.solution.orientation.vel_ii_mag.extend(results[14])

    config.solution.orientation.pos_pp[0].extend(results[15])
    config.solution.orientation.pos_pp[1].extend(results[16])
    config.solution.orientation.pos_pp[2].extend(results[17])
    config.solution.orientation.pos_pp_mag.extend(results[18])
    config.solution.orientation.vel_pp[0].extend(results[19])
    config.solution.orientation.vel_pp[1].extend(results[20])
    config.solution.orientation.vel_pp[2].extend(results[21])
    config.solution.orientation.vel_pp_mag.extend(results[22])

    config.solution.orientation.oe[0].extend(results[23])
    config.solution.orientation.oe[1].extend(results[24])
    config.solution.orientation.oe[2].extend(results[25])
    config.solution.orientation.oe[3].extend(results[26])
    config.solution.orientation.oe[4].extend(results[27])
    config.solution.orientation.oe[5].extend(results[28])

    config.solution.orientation.lat.extend(results[29])
    config.solution.orientation.lon.extend(results[30])
    config.solution.orientation.alt.extend(results[31])
    config.solution.orientation.gamma_ii.extend(results[32])
    config.solution.orientation.gamma_pp.extend(results[33])

    config.solution.orientation.h_ii[0].extend(results[34])
    config.solution.orientation.h_ii[1].extend(results[35])
    config.solution.orientation.h_ii[2].extend(results[36])
    config.solution.orientation.h_pp[0].extend(results[37])
    config.solution.orientation.h_pp[1].extend(results[38])
    config.solution.orientation.h_pp[2].extend(results[39])
    config.solution.orientation.h_ii_mag.extend(results[40])
    config.solution.orientation.h_pp_mag.extend(results[41])

    config.solution.orientation.uD[0].extend(results[42])
    config.solution.orientation.uD[1].extend(results[43])
    config.solution.orientation.uD[2].extend(results[44])
    config.solution.orientation.uE[0].extend(results[45])
    config.solution.orientation.uE[1].extend(results[46])
    config.solution.orientation.uE[2].extend(results[47])
    config.solution.orientation.uN[0].extend(results[48])
    config.solution.orientation.uN[1].extend(results[49])
    config.solution.orientation.uN[2].extend(results[50])
    config.solution.orientation.vN.extend(results[51])
    config.solution.orientation.vE.extend(results[52])
    config.solution.orientation.azi_pp.extend(results[53])


    # Physical Properties
    config.solution.physical_properties.rho.extend(results[54])
    config.solution.physical_properties.T.extend(results[55])
    config.solution.physical_properties.p.extend(results[56])
    config.solution.physical_properties.wind[0].extend(results[57])
    config.solution.physical_properties.wind[1].extend(results[58])
    config.solution.physical_properties.wind[2].extend(results[59])
    config.solution.physical_properties.cL.extend(results[60])
    config.solution.physical_properties.cD.extend(results[61])
    config.solution.physical_properties.aoa.extend(results[62])
    config.solution.physical_properties.S.extend(results[63])

    # Performances
    config.solution.performance.mass.extend(results[64])
    config.solution.performance.heat_rate.extend(results[65])
    config.solution.performance.heat_load.extend(results[66])
    config.solution.performance.T_r.extend(results[67])
    config.solution.performance.q.extend(results[68])

    # Forces
    config.solution.forces.gravity_ii[0].extend(results[69])
    config.solution.forces.gravity_ii[1].extend(results[70])
    config.solution.forces.gravity_ii[2].extend(results[71])
    config.solution.forces.drag_pp[0].extend(results[72])
    config.solution.forces.drag_pp[1].extend(results[73])
    config.solution.forces.drag_pp[2].extend(results[74])
    config.solution.forces.drag_ii[0].extend(results[75])
    config.solution.forces.drag_ii[1].extend(results[76])
    config.solution.forces.drag_ii[2].extend(results[77])
    config.solution.forces.lift_pp[0].extend(results[78])
    config.solution.forces.lift_pp[1].extend(results[79])
    config.solution.forces.lift_pp[2].extend(results[80])
    config.solution.forces.lift_ii[0].extend(results[81])
    config.solution.forces.lift_ii[1].extend(results[82])
    config.solution.forces.lift_ii[2].extend(results[83])
    config.solution.forces.force_ii[0].extend(results[84])
    config.solution.forces.force_ii[1].extend(results[85])
    config.solution.forces.force_ii[2].extend(results[86])
    config.solution.forces.energy.extend(results[87])

    # Simulation
    config.solution.simulation.MC_seed.extend(results[88])
    config.solution.simulation.drag_passage.extend(results[89])
    # runtime = datetime.now() - startTime
    # print("save_results time: ",runtime)
    # with open("/home/karnapp/Desktop/AerobrakingFiles/AerobrakingMaterial/saveres.txt","a") as text_file:
        # text_file.write("\nSave_Results Runtime: {0}".format(runtime))
    return time_0

def clean_results():
    config.solution.orientation.time = []
    config.solution.orientation.year = []
    config.solution.orientation.month = []
    config.solution.orientation.day = []
    config.solution.orientation.hour = []
    config.solution.orientation.min = []
    config.solution.orientation.second = []
    config.solution.orientation.numberofpassage = []
    config.solution.orientation.pos_ii[0] = []
    config.solution.orientation.pos_ii[1] = []
    config.solution.orientation.pos_ii[2] = []
    config.solution.orientation.vel_ii[0] = []
    config.solution.orientation.vel_ii[1] = []
    config.solution.orientation.vel_ii[2] = []
    config.solution.orientation.pos_ii_mag = []
    config.solution.orientation.vel_ii_mag = []

    config.solution.orientation.pos_pp[0] = []
    config.solution.orientation.pos_pp[1] = []
    config.solution.orientation.pos_pp[2] = []
    config.solution.orientation.pos_pp_mag = []
    config.solution.orientation.vel_pp[0] = []
    config.solution.orientation.vel_pp[1] = []
    config.solution.orientation.vel_pp[2] = []
    config.solution.orientation.vel_pp_mag = []

    config.solution.orientation.oe[0] = []
    config.solution.orientation.oe[1] = []
    config.solution.orientation.oe[2] = []
    config.solution.orientation.oe[3] = []
    config.solution.orientation.oe[4] = []
    config.solution.orientation.oe[5] = []

    config.solution.orientation.lat = []
    config.solution.orientation.lon = []
    config.solution.orientation.alt = []
    config.solution.orientation.gamma_ii = []
    config.solution.orientation.gamma_pp = []

    config.solution.orientation.h_ii[0] = []
    config.solution.orientation.h_ii[1] = []
    config.solution.orientation.h_ii[2] = []
    config.solution.orientation.h_pp[0] = []
    config.solution.orientation.h_pp[1] = []
    config.solution.orientation.h_pp[2] = []
    config.solution.orientation.h_ii_mag = []
    config.solution.orientation.h_pp_mag = []

    config.solution.orientation.uD[0] = []
    config.solution.orientation.uD[1] = []
    config.solution.orientation.uD[2] = []
    config.solution.orientation.uE[0] = []
    config.solution.orientation.uE[1] = []
    config.solution.orientation.uE[2] = []
    config.solution.orientation.uN[0] = []
    config.solution.orientation.uN[1] = []
    config.solution.orientation.uN[2] = []
    config.solution.orientation.vN = []
    config.solution.orientation.vE = []
    config.solution.orientation.azi_pp = []

    # Physical Properties
    config.solution.physical_properties.rho = []
    config.solution.physical_properties.T = []
    config.solution.physical_properties.p = []
    config.solution.physical_properties.wind[0] = []
    config.solution.physical_properties.wind[1] = []
    config.solution.physical_properties.wind[2] = []
    config.solution.physical_properties.cL = []
    config.solution.physical_properties.cD = []
    config.solution.physical_properties.aoa = []
    config.solution.physical_properties.S = []

    # Performances
    config.solution.performance.mass = []
    config.solution.performance.heat_rate = []
    config.solution.performance.heat_load = []
    # heat load calculated in the complete passage script # I don't have idea of what this meant
    config.solution.performance.T_r = []
    config.solution.performance.q = []

    # Forces
    config.solution.forces.gravity_ii[0] = []
    config.solution.forces.gravity_ii[1] = []
    config.solution.forces.gravity_ii[2] = []
    config.solution.forces.drag_pp[0] = []
    config.solution.forces.drag_pp[1] = []
    config.solution.forces.drag_pp[2] = []
    config.solution.forces.drag_ii[0] = []
    config.solution.forces.drag_ii[1] = []
    config.solution.forces.drag_ii[2] = []
    config.solution.forces.lift_pp[0] = []
    config.solution.forces.lift_pp[1] = []
    config.solution.forces.lift_pp[2] = []
    config.solution.forces.lift_ii[0] = []
    config.solution.forces.lift_ii[1] = []
    config.solution.forces.lift_ii[2] = []
    config.solution.forces.force_ii[0] = []
    config.solution.forces.force_ii[1] = []
    config.solution.forces.force_ii[2] = []
    config.solution.forces.energy = []

    # Simulation
    config.solution.simulation.MC_seed = []
    config.solution.simulation.drag_passage = []

    # Closed form
    config.solution.closed_form.t_cf = []
    config.solution.closed_form.h_cf = []
    config.solution.closed_form.v_cf = []
    config.solution.closed_form.gamma_cf = []
    return