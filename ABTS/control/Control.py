#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

# from config import *
import math
import numpy as np
import config
from utils.closed_form_solution import closed_form
from physical_models.Density_models import density_exp as dm
from scipy import optimize
import scipy

def nocontrol(ip,m, args=0,rho=0, T_p=0, T_w=0, S=0, position=0, current_position=0, t=0,index_ratio = 0,state=0):
    angle_of_attack = m.aerodynamics.aoa
    return angle_of_attack


def control_solarpanels_heatrate(ip,m,index_ratio,state,args=0,t=0, position=0, current_position=0):
    if index_ratio[0] == 1:
        T_p = state[0]
        rho = state[1]
        S = state[2]
        if args is not 0:
            if args.montecarlo == True:
                from physical_models.MonteCarlo_perturbations import monte_carlo_guidance_environment
                rho, T_p, S = monte_carlo_guidance_environment(rho , T_p, S ,args)
        T_w = T_p

        taf  = m.aerodynamics.thermal_accomodation_factor
        R = m.planet.R
        gamma = m.planet.gamma
        max_angleofattack = m.aerodynamics.aoa#math.pi/2
        min_angleofattack = 0.0001
        thermal_limit = m.aerodynamics.heat_rate_limit - 0.00001

        heat_rate_max = heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, max_angleofattack)

        heat_rate_min = heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, min_angleofattack)
        L = (taf * rho * R * T_p) * ((R * T_p / (2.0 * math.pi)) ** 0.5)*1e-4

        def f ( x ):
            fx = L * ((S ** 2 + (gamma) / (gamma - 1) - (gamma + 1) / (2 * (gamma - 1)) * (T_w / T_p)) * (
                        math.exp(-(S * math.sin(x)) ** 2) + (math.pi ** 0.5) * (S * math.sin(x)) * (
                        1 + math.erf(S * math.sin(x)))) - 0.5 * math.exp(-(S * math.sin(x)) ** 2)) - thermal_limit

            return fx  # W/cm^2

        class FalseException(Exception):
            pass

        if(heat_rate_max < thermal_limit):
            angle_of_attack = max_angleofattack
        elif(heat_rate_min > thermal_limit):
            angle_of_attack = 0.00001
        elif (heat_rate_max >= thermal_limit) and (heat_rate_min <= thermal_limit):  # W/cm^2
            # print('time',t,'heat rate min', heat_rate_min, 'thermal limit', thermal_limit)
            x_0 = config.aoa_past
            try:
                angle_of_attack = optimize.newton(f, x_0, fprime=lambda x: L * S * math.cos(x) * ((math.pi ** 0.5) * (
                    S ** 2 + gamma / (gamma - 1) + (gamma + 1) / (2 * (gamma - 1)) * T_w / T_p) * (1 + math.erf(
                    S * math.sin(x))) + S * math.sin(x) * math.exp(-(S * math.sin(x)) ** 2)), tol=1e-5, maxiter=1000)
                # print(angle_of_attack)
                if angle_of_attack < 0 or angle_of_attack > np.pi / 2:
                # print('here in the if')
                    angle_of_attack = optimize.newton(f, 1e-1, fprime=lambda x: L * S * math.cos(x) * ((math.pi ** 0.5) * (
                        S ** 2 + gamma / (gamma - 1) + (gamma + 1) / (2 * (gamma - 1)) * T_w / T_p) * (1 + math.erf(
                        S * math.sin(x))) + S * math.sin(x) * math.exp(-(S * math.sin(x)) ** 2)), tol=1e-5, maxiter=1000)
                if angle_of_attack < 0 or angle_of_attack > np.pi / 2: raise FalseException()
            except:
                try:
                    if abs(heat_rate_max - thermal_limit) < abs(
                            heat_rate_min - thermal_limit):  # Newton method is unable to find a solution since there are multiple ones. We need to provide a good initial guess
                        x_0 = 2 * max_angleofattack / 3.0
                    elif abs(heat_rate_max - thermal_limit) > abs(heat_rate_min - thermal_limit):
                        x_0 = 2 * max_angleofattack / 6.0
                    angle_of_attack = optimize.newton(f, x_0, fprime=lambda x: L * S * math.cos(x) * ((math.pi ** 0.5) * (
                            S ** 2 + gamma / (gamma - 1) + (gamma + 1) / (2 * (gamma - 1)) * T_w / T_p) * (1 + math.erf(
                        S * math.sin(x))) + S * math.sin(x) * math.exp(-(S * math.sin(x)) ** 2)))
                except:
                    print('Check heat rate controller does not CONVERGE!')
                    angle_of_attack = 0.00001


            # except:
            #     try:
            #         print('heat_rate_max ',heat_rate_max ,'heat_rate_min ',heat_rate_min)
            #         if abs(heat_rate_max - thermal_limit) < abs(
            #                 heat_rate_min - thermal_limit):  # Newton method is unable to find a solution since there are multiple ones. We need to provide a good initial guess
            #             x_0 = 2 * max_angleofattack / 3.0
            #         elif abs(heat_rate_max - thermal_limit) > abs(heat_rate_min - thermal_limit):
            #             x_0 = 2 * max_angleofattack / 6.0
            #         angle_of_attack = optimize.newton(f, x_0, fprime=lambda x: L * S * math.cos(x) * ((math.pi ** 0.5) * (
            #                 S ** 2 + gamma / (gamma - 1) + (gamma + 1) / (2 * (gamma - 1)) * T_w / T_p) * (1 + math.erf(
            #             S * math.sin(x))) + S * math.sin(x) * math.exp(-(S * math.sin(x)) ** 2)))
            #     except:
            #         print('here')
            #         angle_of_attack = config.aoa_past
            #print(math.degrees(x_0),math.degrees(angle_of_attack),'heat rate max',heat_rate_max,'heat rate min',heat_rate_min,'thermal limit',thermal_limit)
        else:
            print('Check controller - second check')


        if (angle_of_attack > max_angleofattack) or (angle_of_attack < 0):
            angle_of_attack = 0
            #print('Warning: No solution for Newton Method for angle of attack - change minimum periapsis!')
        return angle_of_attack
    else:
        return config.aoa

def heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, angle): # if you make a change here, remember to change also the thermal models.py heatrate_convective_maxwellian function
    first_term = np.multiply(1e-4 * taf * R * T_p * (np.sqrt(R * T_p / (2 * math.pi))),
                             rho)
    term_a = np.exp(-np.square(np.multiply(S,np.sin(angle)))) + np.multiply(np.sqrt(math.pi) * (np.multiply(S,np.sin(angle))),
                            (1 + scipy.special.erf(np.multiply(S,np.sin(angle)))))
    term_b = np.multiply((S ** 2 + (gamma) / (gamma - 1) - (gamma + 1) / (2 * (gamma - 1)) * (T_w / T_p)),term_a)

    return np.multiply(term_b- 0.5 * np.exp(
                        -(np.square(np.multiply(S,np.sin(angle))))),first_term )     # W/cm^2


def control_solarpanels_heatload(ip,m,args,index_ratio,t=0,position=0,current_position=0,heat_rate_control=False, state = 0):
    # Import Module
    if args.flash2_through_integration == 1:
        if args.heat_load_sol == 0 or args.heat_load_sol == 1:
            from control.heatload_control.time_switch_calcs import switch_calculation_with_integration as switch_calculation
        if args.heat_load_sol == 2 or args.heat_load_sol == 3:
            from control.heatload_control.second_tsw_calcs import second_time_switch_recalc_with_integration as switch_calculation

        from control.heatload_control.second_tsw_calcs import second_time_switch_recalc_with_integration as second_time_switch_recalc
        args.security_mode = 0
    else:
        from control.heatload_control.time_switch_calcs import switch_calculation
        from control.heatload_control.second_tsw_calcs import second_time_switch_recalc
    from control.heatload_control.security_mode import security_mode

    # IF increase too much the 50 sec, the solution becomes a bit unstable because of Mars Gram
    start_reevaluation = False #t>config.time_switch_1: #
    if (config.ascending_phase) or (config.ascending_phase == False and abs(config.time_switch_2-t)<0.2*config.time_switch_2): # in periapsis passed or the time switch is too close to the current time but periapsis not passed, start reevaluation
        start_reevaluation = True
    # print(config.heat_load_past)
    # Calculate time switches
    if (config.evaluate_switch_heat_load == False):
        [config.time_switch_1,config.time_switch_2] = switch_calculation(ip,m,position, args,t,heat_rate_control, reevaluation_mode = 1, current_position=position)
        # print('time switch 1',config.time_switch_1,'time switch 2',config.time_switch_2)
        config.evaluate_switch_heat_load = True

    #Reevaluate of time switch
    elif (config.evaluate_switch_heat_load == True) and start_reevaluation and((t-config.time_switch_1>1 and t-config.timer_revaluation >10 and t-config.time_switch_2<0) or (3<config.time_switch_2-t<50 and t-config.timer_revaluation >3) or (0<config.time_switch_2-t<3 and t-config.timer_revaluation >0.8)) and config.security_mode == False:
        if (t-config.timer_revaluation) > 3:
            reevaluation_mode =1
        else:
            reevaluation_mode = 2
        config.timer_revaluation = t
        if args.second_switch_reevaluation == 1:
            [config.time_switch_1,config.time_switch_2] = second_time_switch_recalc(ip,m,position,args,t, heat_rate_control, reevaluation_mode = reevaluation_mode, current_position=current_position)
            # print(t, 'Reevaluation -- time switch 1', config.time_switch_1, 'time switch 2', config.time_switch_2)
    # Safe mode if almost reached heat_load_limit

    elif (config.heat_load_past)>0.98*m.aerodynamics.heat_load_limit and config.heat_load_past-config.heat_load_ppast < 2 and args.security_mode == 1 and config.security_mode == False and index_ratio[1] == 1:
        [config.time_switch_1 , config.time_switch_2] = security_mode(ip,m,position,args,t,heat_rate_control=False)
        # print(t, 'Security mode -- time switch 1', config.time_switch_1, 'time switch 2', config.time_switch_2,'heat load past',config.heat_load_past)

    if args.heat_load_sol == 0 or args.heat_load_sol == 3:
        if (t > config.time_switch_1) and (t < config.time_switch_2): #correct
            aoa = 0.0
        else:
            aoa = m.aerodynamics.aoa

    elif args.heat_load_sol == 1 or args.heat_load_sol == 2:
        if (t > config.time_switch_1) and (t < config.time_switch_2):
            aoa = m.aerodynamics.aoa
        else:
            aoa = 0.0
    config.heat_load_ppast = config.heat_load_past
    return aoa


def control_solarpanels_openloop(ip,m, args ,index_ratio, state,t=0,position=0,current_position=0,heat_rate_control=True):
    control_solarpanels_heatload(ip,m,args,index_ratio,t=t,position=position, current_position=current_position, heat_rate_control=True,state=0)
    if args.heat_load_sol == 0 or args.heat_load_sol == 3:
        if t >= config.time_switch_1 and t <= config.time_switch_2:
            aoa = 0.0
            # print(t, math.degrees(aoa),config.time_switch_2)
        else:
            aoa = control_solarpanels_heatrate(ip,m,index_ratio, state,args,t)
    elif args.heat_load_sol == 1 or args.heat_load_sol == 2:
        if t >= config.time_switch_1 and t <= config.time_switch_2:
            aoa = control_solarpanels_heatrate(ip,m,index_ratio, state,args,t)
        else:
            aoa = 0.0
    config.aoa_past = aoa
    # print(t, math.degrees(aoa))

    return aoa