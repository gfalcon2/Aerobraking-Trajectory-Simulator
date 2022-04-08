#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

import math
from utils.Odyssey_maneuver_plan import Odyssey_firing_plan
import config
def no_manuever(t0,thrust_mag,delta_v,args,index_phase_aerobraking):
    thrust = 0
    return thrust

def abms(t0,thrust_mag,delta_v,args,index_phase_aerobraking):
    if index_phase_aerobraking==0:
        thrust = thrust_mag
    else:
        thrust = 0.0
    return thrust

def deceleration_drag_passage(t0,thrust_mag,delta_v,args,index_phase_aerobraking):
    if config.drag_state == True:
        thrust = thrust_mag
    else:
        thrust = 0
    return thrust
