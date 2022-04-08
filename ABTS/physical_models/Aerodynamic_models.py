#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import math
from math import cos
from math import sin
import config
from physical_models.MonteCarlo_perturbations import monte_carlo_aerodynamics

def aerodynamicscoefficient_constant(args, T=0, S=0, body=0, aoa=0, montecarlo = 0):

    CL_body = 0.0

    CD_body = 2*(2.2-0.8)/math.pi * args.angle_of_attack+0.8

    if montecarlo == True:
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)

    return CL_body, CD_body


def aerodynamicscoefficient_fM (aoa, body, T, S, args, montecarlo=0):
    alpha = aoa
    sigma = args.reflection_coefficient
    Tw = T

    def pressure(S,aoa,rho_inf,velocity,sigma):
        p = (rho_inf*velocity**2)/(2*S**2)*((((2-sigma)/math.pi**0.5)*S*sin(aoa)+((T/Tw)**0.5)*sigma/2)*math.exp(-(S*sin(aoa))**2)+((2-sigma)*((S*sin(aoa))**2+1/2)+sigma/2*math.pi**0.5*(S*sin(aoa)))*(1+math.erf(S*sin(aoa))))
        return p

    def tau(S,aoa,rho_inf,velocity,sigma):
        t = ((sigma * math.cos(aoa) * rho_inf*velocity**2) / (math.pi ** 0.5 * 2* S)) * (
                math.exp(-(S * sin(aoa)) ** 2) + math.pi ** 0.5 * (S * sin(aoa)) * (
                    1 + math.erf(S * sin(aoa))))
        return t

    def normalcoefficient(S,aoa,sigma):
        CN = 1 / (S ** 2) * ((((2 - sigma) / math.pi ** 0.5) * S * sin(aoa) + sigma / 2.0) * math.exp(
            -(S * sin(aoa)) ** 2) + (
                                         (2.0 - sigma) * ((S * sin(aoa)) ** 2.0 + 0.5) + sigma / 2.0 * math.pi ** 0.5 * (
                                             S * sin(aoa))) * (1.0 + math.erf(S * sin(aoa))))
        return CN

    def axialcoefficient(S,aoa,sigma):
        CA = ((sigma * cos(aoa)) / (math.pi ** 0.5 * S)) * (
                    math.exp(-(S * sin(aoa)) ** 2) + math.pi ** 0.5 * (S * sin(aoa)) * (
                        1.0 + math.erf(S * sin(aoa))))
        return CA


    ## Solar Panels:
    CN_sa = normalcoefficient(S,alpha,sigma)
    CA_sa = axialcoefficient(S,alpha,sigma)
    CL_sa = CN_sa*cos(alpha)-CA_sa*sin(alpha)
    CD_sa = CA_sa*cos(alpha)+CN_sa*sin(alpha)

    ## Spacecraft
    CN_sc = normalcoefficient(S,math.pi*0.5,sigma)
    CA_sc = axialcoefficient(S,math.pi*0.5,sigma)
    CL_sc = CN_sc*cos(math.pi*0.5)-CA_sc*sin(math.pi*0.5)
    CD_sc = CA_sc*cos(math.pi*0.5)+CN_sc*sin(math.pi*0.5)

    CD_body = (CD_sa*body.Area_SA + CD_sc*body.Area_SC)/(body.Area_SA+body.Area_SC)
    CL_body = (CL_sa*body.Area_SA + CL_sc*body.Area_SC)/(body.Area_SA+body.Area_SC)

    if montecarlo == True:
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)

    return CL_body, CD_body

def aerodynamicscoefficient_noballisticflight(aoa, body, args, T=0, S=0, a=0, montecarlo = 0):
    # Newtonian Aerodynamics
    NoseRadius = body.NoseRadius
    BaseRadius = body.BaseRadius
    k = NoseRadius/BaseRadius
    Cp_max = 2
    delta = math.radians(body.delta)

    #from Bonnet Braun wrong
    # first_term_CN = (k**2)*0.5*sin(2*aoa)*(cos(delta))**4
    # second_term_CN = (1-k**2*(cos(delta))**2)*(cos(delta))**2*sin(2*aoa)
    # CN_body = Cp_max*0.5*(first_term_CN+second_term_CN)
    # first_term_CA = (k**2)*0.5*(2*cos(delta)**2*(1-sin(delta)**4)+sin(aoa)**2*cos(delta)**4)
    # second_term_CA = (1-k**2*cos(delta)**2)*(2*sin(delta)**2*cos(aoa)**2+cos(delta)**2*sin(aoa)**2)
    # CA_body = Cp_max*0.5*(first_term_CA+second_term_CA)

    #from Putnam slides
    CA_body = (1-sin(delta)**4)*k**2+(2*sin(delta)**2*cos(aoa)**2+cos(delta)**2*sin(aoa)**2)*(1-k**2*cos(delta)**2)
    CN_body = (1-k**2*cos(delta)**2)*cos(delta**2)*sin(2*aoa)


    CD_body = CA_body*cos(aoa)+CN_body*sin(aoa)-0.15
    CL_body = CN_body*cos(aoa)-CA_body*sin(aoa)

    if montecarlo == True:
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body,args)

    return CL_body, CD_body

