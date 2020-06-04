#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

import math
import config
def no_manuever(OE,t0,thrust_mag,args):
    thrust = 0
    return thrust

def abms(vi,t0,thrust_mag,args):
    if abs(vi-math.pi)<=0.1 and config.firing_on == 0 and config.firing_orbit == 0:
        config.firing_on=1
        config.firing_time = t0
        config.firing_orbit =1
        thrust = thrust_mag
    elif config.firing_on == 0:
        thrust = 0
    elif config.firing_on==1:
        thrust = thrust_mag
        if (t0 - config.firing_time) >= args.firing_time:
            config.firing_on = 0
    return thrust


def deceleration_drag_passage(vi,t0,thrust_mag,args):
    if config.drag_state == True:
        thrust = thrust_mag
    else:
        thrust = 0
    return thrust
