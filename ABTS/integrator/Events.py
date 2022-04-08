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
    if len(solution.t_events[-1]) != 0:
        print('IMPACT!')
        breaker = False
    elif len(solution.t_events[-2]) != 0: # If a switch between apoapsis and periapsis occur, break loop
        breaker = False
        print('PERIAPSIS GREATER THAN APOAPSIS!')
    return breaker