#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""


def event(solution):
# Impact Definition
    breaker = True
    if solution.t_events[-1]:
        print('Impact!')
        breaker = False
    elif solution.t_events[-2]: # If a switch between apoapsis and periapsis occur, break loop
        breaker = False
        print('Periapsis Greater than Apoapsis!')
    return breaker