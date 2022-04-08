#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

from utils.Reference_system import *
import numpy as np
import config
def gravity_const(pos_ii_mag, pos_ii,p, mass =0,vel_ii=0):
    mu = p.mu
    pos_ii_hat = pos_ii / pos_ii_mag
    if config.drag_state == False: # gravity constant only in drag passage, otherwise inversed squared
        gravity_ii_mag = -mu / pos_ii_mag ** 2
    else:
        gravity_ii_mag= -p.g_ref
    g = gravity_ii_mag * (pos_ii_hat)
    return g

def gravity_invsquared (pos_ii_mag, pos_ii, p, mass =0,vel_ii=0):
    mu = p.mu
    pos_ii_hat = pos_ii / pos_ii_mag
    gravity_ii_mag_spherical = -mu/ pos_ii_mag ** 2
    g = gravity_ii_mag_spherical * (pos_ii_hat)
    return g

def gravity_invsquaredandJ2effect (pos_ii_mag, pos_ii, p, mass, vel_ii=0):
    mu = p.mu
    J2 = p.J2

    method = 3

    if method == 1: # http://control.asu.edu/Classes/MAE462/462Lecture13.pdf
        # calculation of J2 in RTN system and transformation from RTN to inertial
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, p)
        a , e , i , OMEGA , omega , vi = OE[0] , OE[1] , OE[2] , OE[3] , OE[4] , OE[5]
        pos_ii_hat = pos_ii / pos_ii_mag
        gravity_ii_mag_spherical = -mu/ pos_ii_mag ** 2
        u = omega+vi #latitude angle:
        J2_rtn = -3.0*p.mu*J2*p.Rp_e**2/pos_ii_mag**4* np.array([0.5-1.5*np.sin(i)**2* np.sin(u)**2,
                                                           np.sin(i) ** 2 * np.sin(u)*np.cos(u),
                                                            np.sin(i)*np.cos(i)*np.sin(u)])

        T_ijk = [[np.cos(OMEGA) * np.cos(u) - np.sin(OMEGA) * np.sin(u) * np.cos(i),
                  np.sin(OMEGA) * np.cos(u) + np.cos(OMEGA) * np.sin(u) * np.cos(i),
                  np.sin(u) * np.sin(i)],
                 [-np.cos(OMEGA) * np.sin(u) - np.sin(OMEGA) * np.cos(u) * np.cos(i),
                  -np.sin(OMEGA) * np.sin(u) + np.cos(OMEGA) * np.cos(u) * np.cos(i),
                  np.cos(u) * np.sin(i)],
                 [np.sin(OMEGA) * np.sin(i), -np.cos(OMEGA) * np.sin(i), np.cos(i)]]
        J2_ii = np.inner(np.transpose(T_ijk),J2_rtn)
        g = gravity_ii_mag_spherical * (pos_ii_hat) + J2_ii

    elif method == 2: # https://link.springer.com/content/pdf/10.1007/978-981-10-2383-5_2.pdf
        # calculation of J2 in LVLH system and transformation from LVLH to inertial
        OE = rvtoorbitalelement(pos_ii, vel_ii, mass, p)
        a , e , i , OMEGA , omega , vi = OE[0] , OE[1] , OE[2] , OE[3] , OE[4] , OE[5]
        u = omega+vi #latitude angle:
        T_ijk = [[np.cos(OMEGA) * np.cos(u) - np.sin(OMEGA) * np.sin(u) * np.cos(i),
                  np.sin(OMEGA) * np.cos(u) + np.cos(OMEGA) * np.sin(u) * np.cos(i),
                  np.sin(u) * np.sin(i)],
                 [-np.cos(OMEGA) * np.sin(u) - np.sin(OMEGA) * np.cos(u) * np.cos(i),
                  -np.sin(OMEGA) * np.sin(u) + np.cos(OMEGA) * np.cos(u) * np.cos(i),
                  np.cos(u) * np.sin(i)],
                 [np.sin(OMEGA) * np.sin(i), -np.cos(OMEGA) * np.sin(i), np.cos(i)]]

        # transformation of r from ECI to LVLH
        r_LHLV = np.inner(np.transpose(T_ijk),pos_ii)
        r_LHLV_mag = np.linalg.norm(r_LHLV)
        kJ2 = 1.5 *J2*mu*p.Rp_e**2
        grad_U = [mu/r_LHLV_mag**2 + (kJ2/r_LHLV_mag**4)*(1.0-3.0*np.sin(i)**2*np.sin(u)**2),
                  (kJ2/r_LHLV_mag**4)*np.sin(i)**2*np.sin(2.0*u),
                  (kJ2 / r_LHLV_mag ** 4) * np.sin(2.0*i) * np.sin(u)]
        g_LHLV = [-item for item in grad_U]
        # transformation of gravitational acceleration from LVLH to ECI
        g = np.inner(T_ijk,g_LHLV)

    elif method == 3: # https://www.vcalc.com/equation/?uuid=1e5aa6ea-95a3-11e7-9770-bc764e2038f2
        # Calculation of J2 in ECI system
        pos_ii_hat = pos_ii / pos_ii_mag
        gravity_ii_mag_spherical = -mu/ (pos_ii_mag ** 2)
        x = pos_ii[0]
        y = pos_ii[1]
        z = pos_ii[2]
        r = pos_ii_mag
        J2_ii =-1.5*J2*mu*(p.Rp_e**2)/(r**4)* np.array([x/r*(5*(z**2/r**2)-1),
                                                         y/r*(5*(z**2/r**2)-1),
                                                         z/r*(5*(z**2/r**2)-3)])
        # g = gravity_ii_mag_spherical * (pos_ii_hat) + J2_ii
        gx = -mu*x/ (pos_ii_mag ** 3)*(1+3/2*J2*(p.Rp_e/r)**2*(1-5*(z/r)**2))
        gy = -mu*y/ (pos_ii_mag ** 3)*(1+3/2*J2*(p.Rp_e/r)**2*(1-5*(z/r)**2))
        gz = -mu*z/ (pos_ii_mag ** 3)*(1+3/2*J2*(p.Rp_e/r)**2*(3-5*(z/r)**2))
        g = np.array([gx,gy,gz])
    return g

# I tried calculation of J2 in different ways, all of them provide always the same results. From one passage to the other, periapsis altitude doesn't change that much, because since omega and OMEGA are set to be 0, the periapsis always occurs to the equator.
# However, the J2 effects the omega and OMEGA rates.