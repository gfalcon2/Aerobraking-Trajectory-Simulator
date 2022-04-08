#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import numpy as np
from utils.Odyssey_maneuver_plan import Odyssey_firing_plan
import config
from utils.Reference_system import *

def RK4(f,h,t,y,m,T_ijk, index_phase_aerobraking, args, solution):
    y = np.array(y)
    k_1 = f(t , y, m, index_phase_aerobraking)
    k_2 = f(t+h/2, y+k_1*h/2, m, index_phase_aerobraking)
    k_3 = f(t+h/2, y+k_2*h/2, m, index_phase_aerobraking)
    k_4 = f(t+h, y+h*k_3, m, index_phase_aerobraking)

    k_1 = k_1.tolist()
    k_2 = k_2.tolist()
    k_3 = k_3.tolist()
    k_4 = k_4.tolist()

    summation = []
    for i in range(len(k_1)):
        summation.append(k_1[i] + 2*k_2[i] + 2*k_3[i] + k_4[i])

    y_n = [0,0,0,0,0,0,0,0]
    for i in range(len(k_1)):
        y_n[i] = y[i] + h*summation[i]/6
    t_n = t + h

    # Check terminal events
    breaker, impact_breaker, solution = events(t_n, y_n, t, y, m, T_ijk, index_phase_aerobraking, args, solution)
    return y_n,t_n,breaker,solution


def events(t_n, y_n, t_0, y_0, m, T_ijk, index_phase_aerobraking,args, solution):

    while True:
        impact_breaker = False
        #First Event
        breaker, solution= impact(t_n, y_n, m, solution,args)  # Stop simulation if impact occurs

        if breaker == True:
            impact_breaker = True
            break

        #Second Event
        if index_phase_aerobraking == 0:   # 1 step simulation
            breaker = apoasispoint(t_n,y_0, y_n,m,args)

            stop_firing(t_n,y_n,m,args)

        elif index_phase_aerobraking == 2:  # 3 step simulation - drag passage
            if args.control_mode != 0 :
                breaker = heat_check(t_n, y_n, t_0, y_0, m,args)
            elif args.drag_passage == True or args.body_shape == 'Blunted Cone': #if the dragpassage key has been defined and is set to true, stop the simulation when altitude = 160 km
                breaker = out_drag_passage(t_n, y_n, t_0, y_0, m,args)
            else:
                breaker = eventsecondstep(t_n, y_n, t_0, y_0, m,args)


        elif index_phase_aerobraking == 1: # 3 step simulation - Initial outer-atmosphere phase, where the abms occur
            stop_firing(t_n , y_n , m , args)
        if breaker == True:
            break

        # Third Event
        breaker, solution = apoapsisgreaterperiapsis(t_n, y_n, m, solution,args)

        break
    return breaker, impact_breaker, solution

def impact(t,y, m, solution,args):#nota: Event impact has to be the last because of the impact condition in the aerobraking script.
    if ((y[0]**2 + y[1]**2 + y[2]**2)**0.5 <= (m.planet.Rp_e+35)):
        breaker = True
        print('Attention: Impact!!')
        solution.t_events[-1] = ['True']
    else:
        breaker = False
    return breaker, solution


def apoasispoint(t,y_0,y_n,m,args):
    pos_ii = [y_n[0] , y_n[1] , y_n[2]]  # Inertial position
    vel_ii = [y_n[3] , y_n[4] , y_n[5]]  # Inertial velocity
    r_i = cartesian(pos_ii[0] , pos_ii[1] , pos_ii[2])
    v_i = cartesian(vel_ii[0] , vel_ii[1] , vel_ii[2])
    [OE_n] = rvtoorbitalelement(r_i , v_i , y[6] , m.planet)

    pos_ii = [y_0[0] , y_0[1] , y_0[2]]  # Inertial position
    vel_ii = [y_0[3] , y_0[4] , y_0[5]]  # Inertial velocity
    r_i = cartesian(pos_ii[0] , pos_ii[1] , pos_ii[2])
    v_i = cartesian(vel_ii[0] , vel_ii[1] , vel_ii[2])
    [OE_0] = rvtoorbitalelement(r_i , v_i , y[6] , m.planet)
    if (OE_n.vi - math.pi)>0 and (OE_0.vi - math.pi)<0:
        breaker = True
    else:
        breaker = False
    return breaker

def apoapsisgreaterperiapsis(t, y, m,solution,args):
    r = y[0:3]
    v = y[3:6]

    Energy = (np.inner(v, v)) / 2 - m.planet.mu / (np.inner(r, r)) ** 0.5
    a = - m.planet.mu / (2 * Energy)
    h = np.cross(r, v)
    h += 0.
    e = (1 + (2 * Energy * (np.inner(h, h)) / (m.planet.mu) ** 2)) ** 0.5

    r_a = a*(1+e)
    r_p = a*(1-e)
    if r_a < r_p:
        breaker = True
        solution.t_events[-2] = ['True']
        print('Periapsis greater than apoapsis!')
    else:
        breaker = False
    return breaker, solution

def eventsecondstep(t_n, y_n, t_0, y_0, m,args):
    alt_y_n = ((y_n[0]**2 + y_n[1]**2 + y_n[2]**2)**0.5) - m.planet.Rp_e -250*1e3
    alt_y_0 = ((y_0[0]**2 + y_0[1]**2 + y_0[2]**2)**0.5) - m.planet.Rp_e -250*1e3
    if alt_y_0< 0 and alt_y_n >= 0 :
        breaker = True
    else:
        breaker = False
    return breaker

def heat_check(t_n, y_n, t_0, y_0, m,args):
    alt_y_n = ((y_n[0]**2 + y_n[1]**2 + y_n[2]**2)**0.5) - m.planet.Rp_e -120*1e3
    alt_y_0 = ((y_0[0]**2 + y_0[1]**2 + y_0[2]**2)**0.5) - m.planet.Rp_e -120*1e3
    if alt_y_0 < 0 and alt_y_n > 0 :
        breaker = True
    else:
        breaker = False
    return breaker

def out_drag_passage(t_n, y_n, t_0, y_0, m,args):
    alt_y_n = ((y_n[0]**2 + y_n[1]**2 + y_n[2]**2)**0.5) - m.planet.Rp_e -args.AE*1e3
    alt_y_0 = ((y_0[0]**2 + y_0[1]**2 + y_0[2]**2)**0.5) - m.planet.Rp_e -args.AE*1e3
    if alt_y_0 < 0 and alt_y_n > 0 :
        breaker = True
    else:
        breaker = False
    return breaker

def stop_firing(t_n,y_n,m,args):
    mass = y_n[6]
    delta_v = (m.engine.g_e * m.engine.Isp)*np.log(m.body.Mass/mass) # delta-v variation
    delta_v_max = Odyssey_firing_plan(delta_v , args)
    if delta_v >= (delta_v_max + config.delta_v_man):
        config.firing_on = 0
    return delta_v - (delta_v_max+config.delta_v_man)
