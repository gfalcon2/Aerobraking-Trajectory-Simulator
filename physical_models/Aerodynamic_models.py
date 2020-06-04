#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import math
import config
from physical_models.MonteCarlo_perturbations import monte_carlo_aerodynamics

def aerodynamicscoefficient_constant(args, T=0, S=0, body=0, aoa=0, montecarlo = 0):

    CL_body = 0
    CD_body = 2
    if montecarlo == True:
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)

    return CL_body, CD_body


def aerodynamicscoefficient_fM (aoa, body, T, S, args, montecarlo=0):
    alpha = aoa
    a_f = args.accomodation_factor
    Tw = T

    def pressure(S,aoa,rho_inf,velocity,sigma):
        p = (rho_inf*velocity**2)/(2*S**2)*((((2-sigma)/math.pi**0.5)*S*math.sin(aoa)+((T/Tw)**0.5)*sigma/2)*math.exp(-(S*math.sin(aoa))**2)+((2-sigma)*((S*math.sin(aoa))**2+1/2)+sigma/2*math.pi**0.5*(S*math.sin(aoa)))*(1+math.erf(S*math.sin(aoa))))
        return p

    def tau(S,aoa,rho_inf,velocity,sigma):
        t = ((sigma * math.cos(aoa) * rho_inf*velocity**2) / (math.pi ** 0.5 * 2* S)) * (
                math.exp(-(S * math.sin(aoa)) ** 2) + math.pi ** 0.5 * (S * math.sin(aoa)) * (
                    1 + math.erf(S * math.sin(aoa))))
        return t

    def normalcoefficient(S,aoa,sigma):
        CN = 1 / (S ** 2) * ((((2 - sigma) / math.pi ** 0.5) * S * math.sin(aoa) + sigma / 2) * math.exp(
            -(S * math.sin(aoa)) ** 2) + (
                                         (2 - sigma) * ((S * math.sin(aoa)) ** 2 + 1 / 2) + sigma / 2 * math.pi ** 0.5 * (
                                             S * math.sin(aoa))) * (1 + math.erf(S * math.sin(aoa))))
        return CN

    def axialcoefficient(S,aoa,sigma):
        CA = ((sigma * math.cos(aoa)) / (math.pi ** 0.5 * S)) * (
                    math.exp(-(S * math.sin(aoa)) ** 2) + math.pi ** 0.5 * (S * math.sin(aoa)) * (
                        1 + math.erf(S * math.sin(aoa))))
        return CA


    ## Solar Panels:
    CN_sa = normalcoefficient(S,alpha,a_f)
    CA_sa = axialcoefficient(S,alpha,a_f)
    CL_sa = CN_sa*math.cos(alpha)-CA_sa*math.sin(alpha)
    CD_sa = CA_sa*math.cos(alpha)+CN_sa*math.sin(alpha)

    ## Spacecraft
    CN_sc = normalcoefficient(S,math.pi/2,a_f)
    CA_sc = axialcoefficient(S,math.pi/2,a_f)
    CL_sc = CN_sc*math.cos(math.pi/2)-CA_sc*math.sin(math.pi/2)
    CD_sc = CA_sc*math.cos(math.pi/2)+CN_sc*math.sin(math.pi/2)

    CD_body = (CD_sa*body.Area_SA + CD_sc*body.Area_SC)/(body.Area_SA+body.Area_SC)
    CL_body = (CL_sa*body.Area_SA + CL_sc*body.Area_SC)/(body.Area_SA+body.Area_SC)

    if montecarlo == True:
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)

    return CL_body, CD_body

## this is not done correctly
def aerodynamicscoefficient_noballisticflight(aoa, body, args, T=0, S=0, a=0, montecarlo = 0):
    NoseRadius = body.NoseRadius
    BaseRadius = body.BaseRadius
    delta = body.delta
    theta = math.pi / 2 - delta

    # Probe
    CD_nose = 1 - (math.cos(theta))**4
    CL_nose = 0

    CN_conicfrostrum = (1-(NoseRadius/BaseRadius)**2 * (math.cos(delta)**2)) * (math.cos(delta))**2 * math.sin(2*aoa)
    CA_conicfrostrum = (1-(NoseRadius/BaseRadius)**2 * (math.cos(delta)**2)) * (2 * (math.sin(delta)**2) * (math.cos(aoa)**2) + (math.cos(delta)**2)*(math.sin(delta)**2))

    CA_body = (NoseRadius/BaseRadius)**2 * CD_nose + CA_conicfrostrum
    CN_body = (NoseRadius/BaseRadius)**2 * CL_nose + CN_conicfrostrum

    CL_body = CN_body * math.cos(aoa) - CA_body * math.sin(aoa)
    CD_body = CA_body * math.cos(aoa) + CN_body * math.sin(aoa)
    if montecarlo == True:
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body,args)

    return CL_body, CD_body

