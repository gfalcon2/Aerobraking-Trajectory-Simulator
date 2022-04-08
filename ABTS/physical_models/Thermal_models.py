#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import math
from scipy.interpolate import interp1d

# Blunt Body
def heatrate_convective(S,T,m,rho,v,aoa):
    rn = m.body.NoseRadius+0.25
    k = m.planet.k#1.7623*10**(-4)  #kg^0.5/m
    q_conv = (k*math.sqrt(rho/rn)*v**3)*1e-4 #W/cm^2
    return q_conv

def heatrate_radiative(S,T,m,rho,v,aoa):
    rn = m.body.NoseRadius
    C = 4.736*1e4
    b = 1.22
    fv = [2040, 1780, 1550, 1313, 1065, 850, 660, 495, 359, 238, 151, 115, 81, 55, 35, 19.5, 9.7, 4.3, 1.5, 0]
    vf = [16000,15500, 15000, 14500,14000, 13500, 13000, 12500, 12000, 11500, 11000, 10750, 10500, 10250, 10000, 9750, 9500, 9250, 9000, 0] #m/s
    fun = interp1d(vf,fv)
    a = 1.072*1e6*v**(-1.88)*rho**(-0.325)
    f = fun(v)
    q_rad = (C*(rn**a) * (rho**b) * f) #W/cm^2
    return q_rad


def heatrate_convective_radiative(S,T,m,rho,v,aoa):
    q_conv = heatrate_convective(S,T,m,rho,v,aoa)
    q_rad = heatrate_radiative(S,T,m,rho,v,aoa)
    return q_conv#+q_rad

# Flat Plates (from Shaaf and Chambre)
def heatrate_convective_maxwellian(S,T,m,rho,v,aoa):
    # if you make a change here, remember to change also the Control.py heatrate_calcs
    a = m.aerodynamics
    p = m.planet
    alpha = a.aoa
    t_m = a.thermal_contact # for thermal perfect contact between rear and frontal area t_f = 0, for thermal insulation between rear and frontal area t_f = 1,

    # Front and rear surface in perfect thermal contact
    if (t_m != 1):
        r_prime = (1/(S**2)) * (2*S**2 + 1 - 1/(1 + (math.pi**0.5)*S*math.sin(alpha)*math.erf(S*math.sin(alpha)*math.exp((S*math.sin(alpha))**2))))
        St_prime = (1/(4*(math.pi**0.5)*S))*(math.exp(-(S*math.sin(alpha))**2) + (math.pi**0.5)*(S*math.sin(alpha))*math.erf(S*math.sin(alpha)))

    # Front and rear surfaces insulated from one another
    elif (t_m == 1):
        r_prime = (1/(S**2))* (2*S**2 + 1 - 1/(1 + (math.pi**0.5)*S*math.sin(alpha)*(1+math.erf(S*math.sin(alpha)))*math.exp((S*math.sin(alpha))**2)))
        St_prime = (1/(4*(math.pi**0.5)*S))*(math.exp(-(S*math.sin(alpha))**2) + (math.pi**0.5)*(S*math.sin(alpha))*(1+math.erf(S*math.sin(alpha))))

    T_0 = T * (1+ ((p.gamma - 1)/p.gamma)*S**2)
    T_r = T + (p.gamma/(p.gamma +1))* r_prime* (T_0 - T)
    T_p = T
    gamma = p.gamma
    T_w = T_p

    # cp = m.planet.gamma / (m.planet.gamma - 1) * m.planet.R
    # heat_rate = ((m.aerodynamics.thermal_accomodation_factor*(m.planet.gamma + 1) /  m.planet.gamma) * (St*rho * vel_pp_mag * cp * (T_r - T_p)) )* 10**-4 # W/cm^2

    heat_rate = (m.aerodynamics.thermal_accomodation_factor * rho * m.planet.R * T_p) * (
            (m.planet.R * T_p / (2.0 * math.pi)) ** 0.5) * (
                        (S ** 2.0 + (gamma) / (gamma - 1.0) - (gamma + 1.0) / (2.0 * (gamma - 1)) * (T_w / T_p)) * (
                        math.exp(-(S * math.sin(aoa)) ** 2.0) + (math.pi ** 0.5) * (S * math.sin(aoa)) *
                        (1 + math.erf(S * math.sin(aoa)))) - 0.5 * math.exp(-(S * math.sin(aoa)) ** 2.0)) * 1e-4  # W/cm^2

    return heat_rate

