#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import numpy as np

def RK4(f,h,t,y,m,T_ijk, index_phase_aerobraking, args):
    k_1 = h*f(t , y, m, index_phase_aerobraking)
    k_2 = h*f(t+h/2, y+k_1/2, m, index_phase_aerobraking)
    k_3 = h*f(t+h/2, y+k_2/2, m, index_phase_aerobraking)
    k_4 = h*f(t+h, y+k_3, m, index_phase_aerobraking)

    k_1 = k_1.tolist()
    k_2 = k_2.tolist()
    k_3 = k_3.tolist()
    k_4 = k_4.tolist()

    summation = []
    for i in range(len(k_1)):
        summation.append(k_1[i] + 2*k_2[i] + 2*k_3[i] + k_4[i])

    y_n = [0,0,0,0,0,0,0]
    for i in range(len(k_1)):
        y_n[i] = y[i] + summation[i]/6
    t_n = t + h


    # Check terminal events
    breaker, impact_breaker = events(t_n, y_n, t, y, m, T_ijk, index_phase_aerobraking, args)
    return np.array(y_n),t_n,breaker,impact_breaker


def events(t_n, y_n, t_0, y_0, m, T_ijk, index_phase_aerobraking,args):

    while True:
        impact_breaker = True
        #First Event
        breaker = impact(t_n, y_n, m)  # Stop simulation if impact occurs

        if breaker == False:
            impact_breaker = False
            break

        #Second Event
        if index_phase_aerobraking == 0:   # 1 step simulation
            if args.drag_passage == True: #if the dragpassage key has been defined and is set to true, stop the simulation when altitude = 160 km
                breaker = eventsecondstep(t_n, y_n, t_0, y_0, m)
            else:
                breaker = apoasispoint(t_n, y_n, t_0, y_0, T_ijk)
        elif index_phase_aerobraking == 2:  # 3 step simulation - drag passage
            breaker = eventsecondstep(t_n, y_n, t_0, y_0, m)

        if breaker == False:
            break

        # Third Event
        breaker = apoapsisgreaterperiapsis(t_n, y_n, m)

        break
    return breaker, impact_breaker

def impact(t,y, m):#nota: Event impact has to be the last because of the impact condition in the aerobraking script.
    if ((y[0]**2 + y[1]**2 + y[2]**2)**0.5 <= (m.planet.Rp_e+70000)):
        breaker = False
        print('Attention: Impact!!')
    else:
        breaker = True
    return breaker


def apoasispoint(t_n, y_n, t_0, y_0, T_ijk):
    y_pwj_0 = np.inner(T_ijk,y_0[3:6])
    y_pwj_n = np.inner(T_ijk,y_n[3:6])
    if (y_pwj_n[0] >= 0) and (y_pwj_0[0] <= 0) :
        breaker = False
    else:
        breaker = True
    return breaker

def apoapsisgreaterperiapsis(t, y, m):
    r = y[0:3]
    v = y[3:6]

    Energy = (np.inner(v, v)) / 2 - m.planet.mu / (np.inner(r, r)) ** 0.5
    a = - m.planet.mu / (2 * Energy)
    h = np.cross(r, v)
    h += 0.
    e = (1 + (2 * Energy * (np.inner(h, h)) / (m.planet.mu) ** 2)) ** 0.5

    r_a = a*(1+e)
    r_p = a * (1 - e)
    if r_a < r_p:
        breaker = False
        print('Periapsis greater than apoapsis!')
    else:
        breaker = True
    return breaker

def eventsecondstep(t_n, y_n, t_0, y_0, m):
    alt_y_n = ((y_n[0]**2 + y_n[1]**2 + y_n[2]**2)**0.5) - m.planet.Rp_e -160*10**3
    alt_y_0 = ((y_0[0]**2 + y_0[1]**2 + y_0[2]**2)**0.5) - m.planet.Rp_e -160*10**3
    if alt_y_0 < 0 and alt_y_n > 0 :
        breaker = False
    else:
        breaker = True
    return breaker

