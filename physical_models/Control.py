#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

from config import *
import math
import numpy as np
import config
from utils.closed_form_solution import closed_form
from integrator.Integrators import RK4
from physical_models.Density_models import density_exp as dm

def nocontrol(m, args=0,rho=0, T_p=0, T_w=0, S=0, OE=0, t=0):
    angle_of_attack = m.aerodynamics.aoa
    return angle_of_attack


def control_solarpanels_heatrate(m,args=0,rho=0,T_p=0,T_w=0,S=0,t=0, OE=0):
    taf  = m.aerodynamics.thermal_accomodation_factor
    R = m.planet.R
    gamma = m.planet.gamma
    max_angleofattack = math.radians(90)#math.pi/2
    min_angleofattack = 0.02
    thermal_limit = m.aerodynamics.thermal_limit - 0.00001
    noise = 0.0001
    x_0 = math.pi/3

    epsilon = 100
    heat_rate_max = heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, max_angleofattack)
    heat_rate_min = heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, min_angleofattack)
    counter = 0
    if (heat_rate_max >= thermal_limit) and (heat_rate_min <= thermal_limit):  # W/cm^2
        if abs(heat_rate_max - thermal_limit) < abs(heat_rate_min - thermal_limit):  #Newton method is unable to find a solution since there are multiple ones. We need to provide a good initial guess
            x_0 = math.pi/3
        elif abs(heat_rate_max - thermal_limit) > abs(heat_rate_min - thermal_limit):
            x_0 =  math.pi/6
        while (epsilon > 0.0001) and (counter < 500): # if get stuck in a local minimum
            counter += 1
            L = (taf * rho * R* T_p) * ((R * T_p / (2 * math.pi)) ** 0.5) * 10 ** -4
            fx = L * ((S ** 2 + (gamma) / (gamma - 1) - (gamma + 1) / (2 * (gamma - 1)) * (T_w / T_p)) * (math.exp(-(S * math.sin(x_0)) ** 2) + (math.pi ** 0.5) * (S * math.sin(x_0)) * (
                        1 + math.erf(S * math.sin(x_0)))) - 1 / 2 * math.exp(-(S * math.sin(x_0)) ** 2)) - thermal_limit
            fx_prime = L * S * math.cos(x_0) * ((math.pi ** 0.5) * (S ** 2 + gamma / (gamma - 1) + (gamma + 1) / (2 * (gamma - 1)) * T_w / T_p) * (
                    1 + math.erf(S * math.sin(x_0))) + S * math.sin(x_0) * math.exp(-(S * math.sin(x_0)) ** 2))
            x_1 = x_0 - (fx / fx_prime)
            epsilon = abs(x_1 - x_0)
            angle_of_attack = x_1
            x_0 = x_1
            if counter==499:
                print('Check controller - first check')
    elif(heat_rate_max < thermal_limit):
        angle_of_attack = max_angleofattack
    elif(heat_rate_min > thermal_limit):
        angle_of_attack = 0.00001
    else:
        print('Check controller - second check')


    if (angle_of_attack > max_angleofattack) or (angle_of_attack < 0):
        print(math.degrees(angle_of_attack))
        angle_of_attack = 0
        #print('Warning: No solution for Newton Method for angle of attack - change minimum periapsis!')


    # Evaluation Energy - porta a zero il controllo se ra z # ESPEDIENTE THERMAL LOAD
    #if (config.count_numberofpassage != 1) and ((config.solution.orientation.oe[0][-1] * (1 + config.solution.orientation.oe[1][-1])) - terminal_state[1] <= 10000*10**3):
    #    angle_of_attack = 0
    return angle_of_attack

def heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, angle):
    return (taf * rho * R * T_p) * ((R * T_p / (2 * math.pi)) ** 0.5) * ((S ** 2 + (gamma) / (gamma - 1) - (gamma + 1) / (2 * (gamma - 1)) * (T_w / T_p)) * (math.exp(-(S * math.sin(angle)) ** 2) + (math.pi ** 0.5) * (S * math.sin(angle)) *(1 + math.erf(S * math.sin(angle)))) - 1 / 2 * math.exp(-(S * math.sin(angle)) ** 2)) * 10 ** -4

def switch_calculation(m,OE,T,heat_load_lim,args,heat_rate_control=False):
    import scipy
    t_cf, h_cf, gamma_cf, v_cf =closed_form(args,m,OE,T,aoa=math.pi/2,online=True)
    RT = T*m.planet.R
    S = (v_cf/(2*RT)**0.5)
    rho = dm(h=h_cf, p=m.planet)[0]
    def heat_rate_calc(rho,S,t_cf, aoa):
        first_term = np.multiply(10 ** -4 * m.aerodynamics.thermal_accomodation_factor * RT * (np.sqrt(RT / (2 * math.pi))),rho)
        term_a = np.multiply(np.exp(-(S * math.sin(aoa)) ** 2) + (math.pi ** 0.5) * (S * math.sin(aoa)),
                        (1 + scipy.special.erf(S * math.sin(aoa))))
        term_b = np.multiply((S ** 2 + 0.5),term_a)
        heat_rate = np.multiply(term_b- 1 / 2 * np.exp(
                    -(S * math.sin(aoa)) ** 2),first_term )     # W/cm^2
        heat_load = 0
        for index in range(1, len(t_cf)):
            heat_load = heat_load + heat_rate[index] * (t_cf[index] - t_cf[index - 1])
        return heat_load

    heat_load = heat_rate_calc(rho,S,t_cf,(math.pi/2))


    if heat_load <= heat_load_lim:
        switch_time = 0
    else:
        breaker = True
        time = -50
        while breaker:
            time +=50
            index = (np.argmin(np.abs(t_cf-time)))

            if heat_rate_control == False:
                temp = (list(map(heat_rate_calc, [rho[0:index], (rho[index + 1:-1])], [S[0:index], (S[index + 1:-1])],
                                 [t_cf[0:index], t_cf[index + 1:-1]], [0, math.pi / 2])))

                hl = np.sum(temp)
            else:
                half_sim = int(round(len(t_cf)/2))
                index_alt = [(np.argmin(np.abs(h_cf[0:half_sim] - 113*10**3))),half_sim+(np.argmin(np.abs(h_cf[half_sim:len(t_cf)] - 113*10**3)))]
                if time < t_cf[index_alt[0]]:
                    temp = (
                        list(map(heat_rate_calc, [rho[0:index], rho[index + 1:index_alt[0]], rho[index_alt[1]:-1]], [S[0:index], S[index + 1:index_alt[0]], S[index_alt[1]:-1]],
                                 [t_cf[0:index], t_cf[index + 1:index_alt[0]], t_cf[index_alt[1]:-1]], [0, math.pi / 2, math.pi / 2])))
                    hl = sum(temp) + config.heat_rate_limit* (t_cf[index_alt[1]]-t_cf[index_alt[0]])
                elif time>= t_cf[index_alt[0]] and time< t_cf[index_alt[1]]:
                    temp = (
                        list(map(heat_rate_calc, [rho[0:index], rho[index_alt[1]:-1]], [S[0:index], S[index_alt[1]:-1]],
                                 [t_cf[0:index], t_cf[index_alt[1]:-1]], [0, math.pi / 2])))
                    hl = sum(temp) + config.heat_rate_limit* (t_cf[index_alt[1]]-t_cf[index])
                elif time>=t_cf[index_alt[1]]:
                    temp = (
                        list(map(heat_rate_calc, [rho[0:index], rho[index:-1]], [S[0:index], S[index:-1]],
                                 [t_cf[0:index], t_cf[index:-1]], [0, math.pi / 2])))
                    hl = sum(temp)

            if hl<=heat_load_lim:
                for t in range(time-50,time,1):
                    index = (np.argmin(np.abs(t_cf - t)))

                    if heat_rate_control == False:
                        temp = (list(
                            map(heat_rate_calc, [rho[0:index], (rho[index + 1:-1])], [S[0:index], (S[index + 1:-1])],
                                [t_cf[0:index], t_cf[index + 1:-1]], [math.radians(0), math.pi / 2])))
                        hl = np.sum(temp)
                    else:
                        if t < t_cf[index_alt[0]]:
                            temp = (
                                list(map(heat_rate_calc,
                                         [rho[0:index], rho[index + 1:index_alt[0]], rho[index_alt[1]:-1]],
                                         [S[0:index], S[index + 1:index_alt[0]], S[index_alt[1]:-1]],
                                         [t_cf[0:index], t_cf[index + 1:index_alt[0]], t_cf[index_alt[1]:-1]],
                                         [0, math.pi / 2, math.pi / 2])))
                            hl = sum(temp) + config.heat_rate_limit * (t_cf[index_alt[1]] - t_cf[index_alt[0]])
                        elif t >= t_cf[index_alt[0]] and time < t_cf[index_alt[1]]:
                            temp = (
                                list(map(heat_rate_calc, [rho[0:index], rho[index_alt[1]:-1]],
                                         [S[0:index], S[index_alt[1]:-1]],
                                         [t_cf[0:index], t_cf[index_alt[1]:-1]], [0, math.pi / 2])))
                            hl = sum(temp) + config.heat_rate_limit * (t_cf[index_alt[1]] - t_cf[index])
                        elif t >= t_cf[index_alt[1]]:
                            temp = (
                                list(map(heat_rate_calc, [rho[0:index], rho[index:-1]], [S[0:index], S[index:-1]],
                                         [t_cf[0:index], t_cf[index:-1]], [0, math.pi / 2])))
                            hl = sum(temp)

                    if abs(hl-heat_load_lim)<0.3 or hl<heat_load_lim:
                        switch_time = t
                        breaker=False
                        break
    #switch_time = 0
    return switch_time

def control_solarpanels_heatload(m,args,rho=0,T_p=0,T_w=0,S=0,terminal_state=0,t=0,OE=0,heat_rate_control=False):
    heat_load_lim = args.max_heat_load
    if config.evaluate_switch_heat_load == False:
        config.time_switch = switch_calculation(m,OE,T_p,heat_load_lim,args,heat_rate_control)
        config.evaluate_switch_heat_load = True
    if t < config.time_switch:
        aoa = 0
    else:
        aoa = math.pi/2
    return aoa

def control_solarpanels_openloop(m, args ,rho=0,T_p=0,T_w=0,S=0,terminal_state=0,t=0,OE=0,heat_rate_control=True):
    aoa = control_solarpanels_heatload(m,args,t=t,OE=OE,T_p=T_p,heat_rate_control=True)
    if t < config.time_switch:
        aoa = aoa
    else:
        aoa = control_solarpanels_heatrate(m,args,rho,T_p,T_w,S,terminal_state)
    return aoa




