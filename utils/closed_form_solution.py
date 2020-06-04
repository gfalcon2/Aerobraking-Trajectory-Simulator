#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

import config
import math
import numpy as np
from physical_models.Density_models import density_exp as dm
from physical_models.Aerodynamic_models import aerodynamicscoefficient_fM as am

from utils.Reference_system import orbitalelemtorv
def closed_form(args,mission,initialcondition=0,T=0,online=False,aoa=0):
    if online == False:
        if args.type_of_mission == 'Drag Passage':
            step_time = len(config.solution.orientation.time)
            initialcondition = config.model.initialcondition(config.solution.orientation.oe[0][0],config.solution.orientation.oe[1][0],config.solution.orientation.oe[2][0],config.solution.orientation.oe[3][0],config.solution.orientation.oe[4][0],config.solution.orientation.oe[5][0],config.solution.performance.mass[0],0,0,0,0,0,0)
            T = config.solution.physical_properties.T[0]
            aoa = config.solution.physical_properties.aoa[0]
            t0 = config.solution.orientation.time[0]

            t_cf , h_cf , gamma_cf , v_cf = closed_form_calculation(t0,mission, initialcondition, aoa, T, step_time)
            results(t_cf , h_cf , gamma_cf , v_cf)

        elif args.type_of_mission != 'Drag Passage':
            # Calculate the number of orbit
            number_orbits = config.solution.orientation.numberofpassage[-1]
            length_solution = len(config.solution.orientation.numberofpassage)
            t , h , gamma , v = [0]*length_solution,[0]*length_solution,[0]*length_solution,[0]*length_solution


            for i in range(1,int(number_orbits)+1):
                idx_orbit = [idx for idx , val in enumerate(config.solution.orientation.numberofpassage) if val == i]
                alt = [(config.solution.orientation.pos_ii_mag[item] - mission.planet.Rp_e) for item in idx_orbit]
                alt_index = [idx for idx , val in enumerate(alt) if val <= 160*10**3]
                index = alt_index[0]+idx_orbit[0]
                step_time = len(alt_index)
                initialcondition = config.model.initialcondition(config.solution.orientation.oe[0][index] ,
                                                                 config.solution.orientation.oe[1][index] ,
                                                                 config.solution.orientation.oe[2][index] ,
                                                                 config.solution.orientation.oe[3][index] ,
                                                                 config.solution.orientation.oe[4][index] ,
                                                                 config.solution.orientation.oe[5][index] ,
                                                                 config.solution.performance.mass[index] , 0 , 0 , 0 , 0 , 0 ,
                                                                 0)
                T = config.solution.physical_properties.T[index]
                aoa = config.solution.physical_properties.aoa[index]
                t0 = config.solution.orientation.time[index]

                t_cf , h_cf , gamma_cf , v_cf = closed_form_calculation(t0, mission , initialcondition , aoa , T , step_time)

                t[alt_index[0]+idx_orbit[0]:alt_index[-1]+idx_orbit[0]] = t_cf
                h[alt_index[0] + idx_orbit[0]:alt_index[-1] + idx_orbit[0]] = h_cf
                gamma[alt_index[0] + idx_orbit[0]:alt_index[-1] + idx_orbit[0]] = gamma_cf
                v[alt_index[0] + idx_orbit[0]:alt_index[-1] + idx_orbit[0]] = v_cf
            # For loop for the number of orbits
            results(t , h , gamma , v)
    else:
        t_cf , h_cf , gamma_cf , v_cf = closed_form_calculation(0,mission, initialcondition, aoa,T, step_time=500)

    return t_cf , h_cf , gamma_cf , v_cf

    #drag_passage -> save results and initial conditions given
    #all_passage -> save results and initial conditions not given
    #online -> not save results and initial conditions given

def closed_form_calculation(t0, mission, initialcondition, aoa, T,step_time):
    [pos_ii, vel_ii] = orbitalelemtorv(initialcondition, mission.planet)
    pos_ii = [pos_ii.x, pos_ii.y, pos_ii.z]
    vel_ii = [vel_ii.x, vel_ii.y, vel_ii.z]
    r0 = np.linalg.norm(pos_ii)  # Inertial position magnitude
    v0 = np.linalg.norm(vel_ii)  # Inertial velocity magnitude
    h0 = r0 - mission.planet.Rp_e


    h_ii = np.cross(pos_ii, vel_ii)
    arg = np.median([-1, 1, np.linalg.norm(h_ii) / (r0 * v0)])  # limit to[-1, 1]
    gamma0 = math.acos(arg)
    if np.inner(pos_ii, vel_ii) < 0:
        gamma0 = -gamma0


    initial_state_angle = initialcondition.vi
    e = initialcondition.e
    a = initialcondition.a
    final_state_angle = -initial_state_angle
    E_initialstate = 2 * math.atan(((1 - e) / (1 + e)) ** 0.5 * math.tan((initial_state_angle) / 2)) # eccentric anomaly
    E_finalstate  = 2 * math.atan(((1 - e) / (1 + e)) ** 0.5 * math.tan((final_state_angle) / 2)) # eccentric anomaly

    # evaluate time to reach next state
    delta_t = (a ** 3 / mission.planet.mu) ** 0.5 * ((E_finalstate - e * math.sin(E_finalstate)) - (E_initialstate - e * math.sin(E_initialstate)))
    t_p = delta_t/2
    t_cf = np.linspace(0,delta_t,step_time)

    cost_3 = v0*gamma0

    h_cf = h0 +cost_3*(t_cf-(np.square(t_cf)/(2*t_p)))

    rho= dm(h=h_cf, p=mission.planet)[0]

    RT = T*mission.planet.R
    S = (v0/(2*RT)**0.5)
    CL,CD = am (aoa, mission.body, T, S, mission.aerodynamics)
    CD0 = CD/aoa
    CL0 = CL
    Area = mission.body.Area_SC+mission.body.Area_SA
    mass = initialcondition.m
    Rp = mission.planet.Rp_e

    cost_1 = (rho*CD0*Area)/(2*mass)
    cost_2 = (rho*CL0*Area)/(2*mass)

    a0 = 0.0016  # 0.0034
    c0 = 5e-06  # 2.9648e-06
    mean_a = 3.38  # 3.5839
    mean_c = 2.6  # 2.5413
    mean_b = -8.25  # -14.5169
    mean_d = -0.001  # 2.5413
    # a0 = 0.00007#0.0034
    # c0 = 1.9648e-07#2.9648e-06
    # mean_a = 0.346#3.5839
    # mean_c = 0.25#2.5413
    # mean_b = -0.8#-14.5169
    # mean_d = -0.01#2.5413

    f1 = -0.005 * v0 + 27.87  # -v0*gamma0/t_p#-0.0053*v0+27.8665
    f2 = a0 * (mean_a ** (2 * abs(math.degrees(gamma0) + 3)) * math.exp(mean_b * (v0 / 1000 - 3.7))) + c0 * (
                mean_c ** (2 * abs(math.degrees(gamma0) + 3)) * math.exp(mean_d * (v0 / 1000 - 3.7)))
    # f = (a1 *math.exp(b1 * abs(math.degrees(gamma0) + 3)) * math.exp(mean_b * (v0 / 1000 - 3.7))) + a2 * math.exp((b2 * abs(math.degrees(gamma0) + 3)) * math.exp(mean_d * (v0 / 1000 - 3.7)))
    # f = a0*(mean_a**(2*abs(math.degrees(gamma0)+3))*math.exp(mean_b*(v0/1000-3.7))) + c0*(mean_c**(2*abs(math.degrees(gamma0)+3))*math.exp(mean_d*(v0/1000-3.7)))
    epsilon = f1 + f2 * (t_cf) / (2 * t_p) * aoa

    k1 = (cost_2 + np.divide(1 , (Rp + h_cf)))
    k2 = np.multiply(cost_1 * aoa * cost_3 , (1 - t_cf / t_p))
    k3 = -mission.planet.g_ref - epsilon
    cost = v0 - (k2[0] / k1[0] - ((k2[0] / k1[0]) ** 2 - 4 * (k3[0] / k1[0])) ** 0.5) / 2
    v_cf = (np.divide(k2 , k1) - np.sqrt(
        np.square(np.divide(k2 , k1)) - 4 * np.divide(k3 , k1))) / 2 + cost

    gamma_cf = cost_3 * np.divide((1 - t_cf / t_p) , v_cf)
    t_cf = np.array([item + t0 for item in t_cf])
    return t_cf, h_cf,gamma_cf,v_cf

def results(t_cf,h_cf,gamma_cf,v_cf):
     # Save Results
    config.solution.closed_form.t_cf.extend(t_cf)
    config.solution.closed_form.h_cf.extend(h_cf)
    config.solution.closed_form.gamma_cf.extend(gamma_cf)
    config.solution.closed_form.v_cf.extend(v_cf)