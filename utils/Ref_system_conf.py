#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

class OE:
    def __init__(self, a, e, i, OMEGA, omega, vi, mass):
        self.a = a
        self.e = e
        self.i = i
        self.OMEGA = OMEGA
        self.omega = omega
        self.vi = vi
        self.m = mass

class cartesian:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z



class R_RA_DEC:
    def __init__(self, r, RA, dec):
        self.r = r
        self.RA = RA
        self.dec = dec

class H_LAN_LON:
    def __init__(self, h, LAT, LON):
        self.h = h
        self.LAT = LAT
        self.LON = LON

class uDuNuE:
    def __init__(self, uD, uN, uE):
        self.uD = uD
        self.uN = uN
        self.uE = uE

class clock:
    def __init__(self,year,month,day,hour,minute,second):
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second