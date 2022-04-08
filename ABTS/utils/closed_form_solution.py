#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

import config
from datetime import *
import math
import numpy as np
from physical_models.Density_models import marsgram as dm #density_exp
from physical_models.Density_models import density_exp as dm1 #density_exp
from utils.misc import find_neighbour
from physical_models.Aerodynamic_models import aerodynamicscoefficient_fM as am

from utils.Reference_system import *
def closed_form(args,mission,initialcondition=0,T=0,online=False,aoa=0,aoa_profile=[]):
    if args.body_shape == 'Blunted Cone':
        len_sol = len(config.solution.orientation.time)
        results([0] * len_sol , [0] * len_sol , [0] * len_sol , [0] * len_sol)
        return [0]*len_sol,[0]*len_sol,[0]*len_sol,[0]*len_sol

    if online == False: #post process after passage
        if args.type_of_mission == 'Drag Passage':
            step_time = len(config.solution.orientation.time)
            initialcondition = [config.solution.orientation.oe[0][0],config.solution.orientation.oe[1][0],config.solution.orientation.oe[2][0],config.solution.orientation.oe[3][0],config.solution.orientation.oe[4][0],config.solution.orientation.oe[5][0],config.solution.performance.mass[0]]
            T = config.solution.physical_properties.T[0]
            aoa = config.solution.physical_properties.aoa[0]
            t0 = config.solution.orientation.time[0]

            t_cf , h_cf , gamma_cf , v_cf = closed_form_calculation(args,t0,mission, initialcondition, aoa, T, step_time)
            results(t_cf , h_cf , gamma_cf , v_cf)

        elif args.type_of_mission != 'Drag Passage':
            # Calculate the number of orbit
            number_orbits = config.solution.orientation.numberofpassage[-1]
            length_solution = len(config.solution.orientation.numberofpassage)
            t , h , gamma , v = [0]*length_solution,[0]*length_solution,[0]*length_solution,[0]*length_solution


            for i in range(1,int(number_orbits)+1):
                idx_orbit = [idx for idx , val in enumerate(config.solution.orientation.numberofpassage) if val == i]
                alt = [(config.solution.orientation.pos_ii_mag[item] - mission.planet.Rp_e) for item in idx_orbit]
                alt_index = [idx for idx , val in enumerate(alt) if val <= 160*1e3]
                if len(alt_index) == 0:
                    len_sol = len(config.solution.orientation.time)
                    results([0]*len_sol,[0]*len_sol,[0]*len_sol,[0]*len_sol)
                    return [0]*len_sol,[0]*len_sol,[0]*len_sol,[0]*len_sol
                index = alt_index[0]+idx_orbit[0]
                step_time = len(alt_index)
                initialcondition = [config.solution.orientation.oe[0][index] ,
                                                                 config.solution.orientation.oe[1][index] ,
                                                                 config.solution.orientation.oe[2][index] ,
                                                                 config.solution.orientation.oe[3][index] ,
                                                                 config.solution.orientation.oe[4][index] ,
                                                                 config.solution.orientation.oe[5][index] ,
                                                                 config.solution.performance.mass[index] ]
                T = config.solution.physical_properties.T[index]
                aoa = config.solution.physical_properties.aoa[index]
                t0 = config.solution.orientation.time[index]

                t_cf , h_cf , gamma_cf , v_cf = closed_form_calculation(args,t0, mission , initialcondition , aoa , T , step_time=step_time)

                t[alt_index[0]+idx_orbit[0]:alt_index[0]+len(alt_index)+idx_orbit[0]] = t_cf
                h[alt_index[0] + idx_orbit[0]:alt_index[0]+len(alt_index) + idx_orbit[0]] = h_cf
                gamma[alt_index[0] + idx_orbit[0]:alt_index[0]+len(alt_index) + idx_orbit[0]] = gamma_cf
                v[alt_index[0] + idx_orbit[0]:alt_index[0]+len(alt_index) + idx_orbit[0]] = v_cf
            # For loop for the number of orbits
            results(t , h , gamma , v)
    else: #online for control
        if args.montecarlo == True and config.closed_form_solution_off:
            config.closed_form_solution_off = 0
            state = {'ra':0,'rp':0,'i':0,'OMEGA':0,'omega':0,'vi':0}
            state['ra'] , state['rp'] , state['i'] , state['OMEGA'] , state['omega'] , state['vi'] = initialcondition[0]*(1+initialcondition[1]) , initialcondition[0]*(1-initialcondition[1]) , initialcondition[2] , initialcondition[3] , initialcondition[4] , initialcondition[5]
            from physical_models.MonteCarlo_perturbations import monte_carlo_guidance_closedform
            state = monte_carlo_guidance_closedform(state,args)
            initialcondition[0], initialcondition[1], initialcondition[2], initialcondition[3] , initialcondition[4] , initialcondition[5] = (state['ra']+state['rp'])/2,(state['ra']-state['rp'])/(state['ra']+state['rp']),state['i'] , state['OMEGA'] , state['omega'] , state['vi']
        t_cf , h_cf , gamma_cf , v_cf = closed_form_calculation(args,0,mission, initialcondition, aoa,T,aoa_profile=aoa_profile, step_time=0)
    return t_cf , h_cf , gamma_cf , v_cf

    #drag_passage -> save results and initial conditions given
    #all_passage -> save results and initial conditions not given
    #online -> not save results and initial conditions given

def closed_form_calculation(args,t0, mission, initialcondition, aoa, T,step_time=0,aoa_profile=[],online=0):
    if config.count_numberofpassage != 1:
        t_prev = config.solution.orientation.time[-1]
    else:
        t_prev = mission.initialcondition.time_rot
    [pos_ii_org, vel_ii_org] = orbitalelemtorv(initialcondition, mission.planet)
    pos_ii = pos_ii_org
    vel_ii = vel_ii_org
    r0 = np.linalg.norm(pos_ii)  # Inertial position magnitude
    v0 = np.linalg.norm(vel_ii)  # Inertial velocity magnitude
    h0 = r0 - mission.planet.Rp_e
    # print('h_in',h0)
    [pos_pp, vel_pp] = r_intor_p(pos_ii, vel_ii, mission.planet, t0,
                                 t_prev)  # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
    LatLong = rtolatlong(pos_pp, mission.planet)
    lat = LatLong[1]
    lon = LatLong[2]
    h0 = LatLong[0]
    # print('alt', h0)

    h_ii = np.cross(pos_ii, vel_ii)
    arg = np.median([-1, 1, np.linalg.norm(h_ii) / (r0 * v0)])  # limit to[-1, 1]
    gamma0 = math.acos(arg)
    if np.inner(pos_ii, vel_ii) < 0:
        gamma0 = -gamma0


    initial_state_angle = initialcondition[5]
    e = initialcondition[1]
    a = initialcondition[0]
    final_state_angle = -initial_state_angle
    E_initialstate = 2 * math.atan(((1 - e) / (1 + e)) ** 0.5 * math.tan((initial_state_angle) / 2)) # eccentric anomaly
    E_finalstate  = 2 * math.atan(((1 - e) / (1 + e)) ** 0.5 * math.tan((final_state_angle) / 2)) # eccentric anomaly



    # evaluate time to reach next state
    delta_t = (a ** 3 / mission.planet.mu) ** 0.5 * ((E_finalstate - e * math.sin(E_finalstate)) - (E_initialstate - e * math.sin(E_initialstate)))
    t_p = delta_t/2
    if h0<160*1e3: #if initial condition are lower than drag passage initial condition #this happens only running MC cases
        # let's calculate pos_ii,v_ii for the point of trajectory corresponding to h = 160 km
        h0 = 160*1e3
        r = mission.planet.Rp_e + h0
        OE = rvtoorbitalelement(pos_ii_org,vel_ii_org,mission,mission.planet)
        a, e, i, OMEGA, omega, vi = OE[0], OE[1], OE[2], OE[3], OE[4], OE[5]
        vi = 2*math.pi - math.acos(((a*(1-e**2)/r)-1)/e)
        E_real_finalstate = 2 * math.atan(
            ((1 - e) / (1 + e)) ** 0.5 * math.tan((-vi) / 2))  # eccentric anomaly
        delta_t = (a ** 3 / mission.planet.mu) ** 0.5 * ((E_real_finalstate - e * math.sin(E_real_finalstate)) - (E_initialstate - e * math.sin(E_initialstate)))
        t_p = delta_t/2

    if step_time == 0:
        temp = delta_t*args.trajectory_rate/10
        if temp > len(config.heat_rate_list):
            step_time = temp
        else:
            step_time = len(config.heat_rate_list)
    t_cf = np.linspace(0,delta_t,int(step_time))

    cost_3 = v0*gamma0

    h_cf = h0 +cost_3*(t_cf-(np.square(t_cf)/(2*t_p)))
    rho= dm1(h=h_cf, p=mission.planet)[0]
    # print(rho)
    # date_initial = datetime(year=mission.initialcondition.year, month=mission.initialcondition.month, day=mission.initialcondition.day,
    #                         hour=mission.initialcondition.hour, minute=mission.initialcondition.min,
    #                         second=mission.initialcondition.second)
    # OE = rvtoorbitalelement(pos_ii, vel_ii, mission.body.Mass, mission.planet)
    # time_real = date_initial + timedelta(seconds=t0)
    # timereal = clock(time_real.year, time_real.month, time_real.day, time_real.hour, time_real.minute, time_real.second)
    # rho= dm(h =h_cf, p = mission.planet, OE = OE, lat = lat, lon = lon, timereal = timereal, t0 = t0, tf_prev = t_prev,
    # montecarlo = 0, Wind = 0, args = args)
    # print(rho)

    RT = T*mission.planet.R
    S = (v0/(2*RT)**0.5)
    CL90,CD90 = am (math.pi/2, mission.body, T, S, mission.aerodynamics)
    CL0,CD0 = am (0, mission.body, T, S, mission.aerodynamics)
    Area_tot = mission.body.Area_SC+mission.body.Area_SA
    mass = initialcondition[-1]
    Rp = mission.planet.Rp_e

    if not aoa_profile:
        aoa_profile = [aoa]*len(t_cf)
    else:
        if len(aoa_profile)>len(t_cf):
            aoa_profile = aoa_profile[0:len(t_cf)]
        elif len(aoa_profile)<len(t_cf):
            last_aoa = aoa_profile[-1]
            aoa_profile = aoa_profile+[last_aoa]*(len(t_cf)-len(aoa_profile))

    CD_t = np.add(CD0, np.divide(np.multiply(aoa_profile,(CD90-CD0)),math.pi/2))
    CL_t = np.add(CL0, np.divide(np.multiply(aoa_profile,(CL90-CL0)),math.pi/2))
    cost_1 = np.divide(rho*CD_t*Area_tot,2*mass)
    cost_2 = np.divide(rho*CL_t*Area_tot,2*mass)

    a0 = 0.0016  # 0.0034
    c0 = 5e-06  # 2.9648e-06
    mean_a = 3.38  # 3.5839
    mean_c = 2.6  # 2.5413
    mean_b = -8.25  # -14.5169
    mean_d = -0.001  # 2.5413


    f1 = -0.005 * v0 + 27.87
    f2 = (a0 * (mean_a ** (2 * abs(math.degrees(gamma0) + 3)) * math.exp(mean_b * (v0 / 1000 - 3.7))) + c0 * (
                mean_c ** (2 * abs(math.degrees(gamma0) + 3)) * math.exp(mean_d * (v0 / 1000 - 3.7))))* (t_cf) / (2 * t_p)

    f2_solar_panels = f2 * aoa*mission.body.Area_SA/ Area_tot
    f2_spacecraft = f2 * math.pi/2*mission.body.Area_SC/ Area_tot

    epsilon = f1 + f2_solar_panels+f2_spacecraft

    k1 = (cost_2 + np.divide(1 , (Rp + h_cf)))
    k2 = np.multiply(cost_1 * cost_3 , (1 - t_cf / t_p))
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