#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

from math import cos as c
from math import sin as s
from math import acos as ac
from math import asin
from math import atan as at
from math import atan2
from utils.Ref_system_conf import *
import numpy as np

def r_intor_p(r_i, v_i, planet, t,t_prev):
    # From PCI (planet centered inertial) to PCPF (planet centered/planet fixed)
    rot_angle = np.linalg.norm(planet.omega)*(t+t_prev) # rad
    # print(t,rot_angle)
    L_pi = [[ c(rot_angle), s(rot_angle), 0],
        [-s(rot_angle) , c(rot_angle), 0],
        [0,             0,             1]]

    r_p = np.inner(L_pi, r_i)

    #V_pf_pci = v_i - np.cross(planet.omega, r_i);  #planet relative velocity in PCI frame
    v_p = np.inner(L_pi,(v_i - np.cross(planet.omega, r_i)))

    return r_p, v_p

def r_pintor_i(r_p, v_p, planet, t, t_prev):
    # From PCPF (planet centered/planet fixed) to PCI (planet centered inertial)
    rot_angle = -np.linalg.norm(planet.omega)*(t+t_prev) # rad
    L_pi = [[ c(rot_angle), s(rot_angle), 0],
        [-s(rot_angle) , c(rot_angle), 0],
        [0,             0,             1]]

    r_i = np.inner((L_pi), r_p)
    v_i = np.inner((L_pi),(v_p + np.cross(planet.omega, r_p)))
    return r_I, v_I



def orbitalelemtorv(oe,planet):
    # From orbital element to ECI (Planet Centered Inertial)

    a,e,i,OMEGA,omega,vi = oe[0],oe[1],oe[2],oe[3],oe[4],oe[5]
    p = a*(1-e**2)
    h = (planet.mu*p)**0.5
    r_x = (h**2)/planet.mu*(1/(1+e*c(vi)))*np.array([c(vi),s(vi),0])
    v_x = planet.mu/h*np.array([-s(vi),e+c(vi),0])
    Q = [[-s(OMEGA)*c(i)*s(omega)+c(OMEGA)*c(omega), c(OMEGA)*c(i)*s(omega)+s(OMEGA)*c(omega), s(i)*s(omega)],
         [-s(OMEGA)*c(i)*c(omega)-c(OMEGA)*s(omega), c(OMEGA)*c(i)*c(omega)-s(OMEGA)*s(omega), s(i)*c(omega)],
         [s(OMEGA)*s(i),                             -c(OMEGA)*s(i),                            c(i)]]

    r = np.matmul(np.transpose(Q),r_x)
    v = np.matmul(np.transpose(Q),v_x)


    # matotor = [c(vi)*(c(OMEGA)*c(omega) - s(OMEGA)*s(omega)*c(i)) - s(vi)*(c(OMEGA)*s(omega) + s(OMEGA)*c(omega)*c(i)),
    #            c(vi)*(s(OMEGA)*c(omega) + c(OMEGA)*s(omega)*c(i)) - s(vi)*(s(OMEGA)*s(omega) - c(OMEGA)*c(omega)*c(i)),
    #            c(vi)*s(omega)*s(i) + s(vi)*c(omega)*s(i)]
    # [x,y,z] = map(lambda x: (p/(1+e*c(vi)))*x, matotor)
    #
    # matotov = [-s(vi)*(c(OMEGA)*c(omega) - s(OMEGA)*s(omega)*c(i)) - (e + c(vi))*(c(OMEGA)*s(omega) + s(OMEGA)*c(omega)*c(i)),
    #            -s(vi)*(s(OMEGA)*c(omega) + c(OMEGA)*s(omega)*c(i)) - (e + c(vi))*(s(OMEGA)*s(omega) - c(OMEGA)*c(omega)*c(i)),
    #            -s(vi)*s(omega)*s(i) + (e + c(vi))*c(omega)*s(i)]
    # matotov = [x+0. for x in matotov]
    # [vx, vy, vz] = map(lambda x: x * (planet.mu/p)**(0.5), matotov)
    # R = [x,y,z]
    # V = [vx, vy, vz]
    R = r.tolist()
    V = v.tolist()
    return R, V

def rvtoorbitalelement(r,v,m,planet):
    # From PCI (Planet Centered Inertial) to orbital element
    i_x = [1,0,0]
    i_y = [0,1,0]
    i_z = [0,0,1]


    Energy = (np.inner(v, v))/2 - planet.mu/(np.inner(r, r))**0.5
    a = - planet.mu/(2* Energy)
    h = np.cross(r, v)
    h += 0.
    index = 0

    r_ver = r/np.linalg.norm(r)
    e = np.cross(v,h)/planet.mu - r_ver
    #e = (1+ (2*Energy*(np.inner(h,h))/(planet.mu)**2))**0.5
    i = ac(np.inner(i_z,h)/np.linalg.norm(h))
    e_vers = e / np.linalg.norm(e)
    if (i == 0) or (i == np.pi):
        if np.inner(e,i_y) >= 0:
            periapsis_longitude = ac(np.inner(i_x,e_vers))
        elif np.inner(e,i_y) < 0:
            periapsis_longitude = 2*np.pi - ac(np.inner(i_x,e_vers))
        OMEGA = periapsis_longitude
        omega = 0

    else:
        n = np.cross(i_z,h)/np.linalg.norm(np.cross(i_z,h))
        if np.inner(n,i_y) >= 0:
            OMEGA = ac(np.inner(i_x,n))
        elif np.inner(n,i_y) < 0:
            OMEGA = 2*np.pi - ac(np.inner(i_x,n))

        if np.inner(e,i_z) >= 0 :
            if (np.inner(n,e_vers) > 1 and np.inner(n,e_vers) < 1+1e-4):
                omega = ac(1)
            elif (np.inner(n,e_vers) < -1 and np.inner(n,e_vers) > -(1+1e-4)):
                omega = ac(-1)
            else:
                omega = ac(np.inner(n,e_vers))
        elif np.inner(e,i_z) < 0:
            if np.inner(n,e_vers) > 1 and np.inner(n,e_vers) < 1+1e-4:
                omega = 2*np.pi - ac(1)
            elif (np.inner(n, e_vers) < -1 and np.inner(n, e_vers) > -(1 + 1e-4)):
                omega = 2*np.pi - ac(-1)
            else:
                omega = 2 * np.pi - ac(np.inner(n, e_vers))

    if np.inner(r,v) > 0:
        value = np.inner(e_vers,r_ver)
        if abs(value) >= 1:
            value = round(value)
        vi = ac(value)
    elif np.inner(r,v) <= 0:
        value = np.inner(e_vers,r_ver)
        if abs(value) >= 1:
            value = round(value)
        vi = 2*np.pi - ac(value)

    e = np.linalg.norm(e)
    if OMEGA == np.pi:
        OMEGA = 0
    return [a,e,i,OMEGA,omega,vi,m]

def rtoalfadeltar(r):
    # From PCI (Planet Centered Inertial) to Geocentric Celestial Reference Frame (GCRF)
    #Conversion between x,y,z and right ascension (RA), declination (dec), distance from the center of the planet (r)
    x = r[0]
    y = r[1]
    z = r[2]
    r = (x**2 + y**2 +z**2)**(0.5)
    l,m,n = x/r, y/r, z/r

    dec = asin(n)
    if m > 0:
        RA = ac(l/c(dec))
    else:
        RA = 2*np.pi - ac(l/c(dec))

    return [r,RA,dec]


def alfadeltartor(R_RA_DEC):
    # From Geocentric Celestial Reference Frame (GCRF) to  PCI (Planet Centered Inertial)
    R = R_RA_DEC[0]
    RA = R_RA_DEC[1]
    DEC = R_RA_DEC[2]

    x = R * c(DEC) *c(RA)
    y = R * c(DEC) *s(RA)
    z = R * s(DEC)
    return [x,y,z]

def latlongtor(LATLONGH, planet, alfa_g0, t, t0): #maybe need to add t_prev orbit
    #From Geodetic to PCI
    phi = LATLONGH[0]  #geodetic latitude
    lam = LATLONGH[1] #longitude
    h = LATLONGH[2] #altitude

    a = planet.Rp_e
    b = planet.Rp_p
    e = (1- b**2/a**2)**0.5
    alfa = lam + alfa_g0 + planet.omega[2]*(t-t0)
    const = a/(1-e**2*s(phi)**2) + h
    x = const * c(phi) * c(alfa)
    y = const * c(phi) * s(alfa)
    z = const * s(phi)

    return [x,y,z]


def rtolatlong(r_p,planet):
    # From PCPF to LLA through Bowring's method https://www.mathworks.com/help/aeroblks/ecefpositiontolla.html;jsessionid=2ae36964c7d5f2115d2c21286db0?nocookie=true
    x_p = r_p[0]
    y_p = r_p[1]
    z_p = r_p[2]

    f = (planet.Rp_e-planet.Rp_p)/planet.Rp_e # flattening
    e = (1-(1-f)**2) # ellipticity (NOTE =  considered as square)
    r = (x_p**2 + y_p**2)**0.5


    # Calculate initial guesses for reduced latitude (latr) and planet-detic latitude (latd)
    latr = at(z_p / ( (1-f)*r ))  # reduced latitude
    latd = at((z_p + (e*(1-f)*planet.Rp_e*s(latr)**3)/(1-e)) / ( r - e*planet.Rp_e*c(latr)**3 ) )

    # Recalculate reduced latitude based on planet-detic latitude
    latr2 = at( (1-f)*s(latd) / c(latd) )
    diff = latr - latr2

    # Iterate until reduced latitude converges
    while diff > 1e-10:
        latr = latr2
        latd = at((z_p + (e*(1-f)*planet.Rp_e*s(latr)**3)/(1-e)) / ( r - e*planet.Rp_e*c(latr)**3 ) )
        latr2 = at( (1-f)*s(latd) /  c(latd) )
        diff = latr - latr2
    lat = latd

    #Calculate longitude
    lon = atan2(y_p,x_p) # -180<lon<180
    #Calculate altitude
    N = planet.Rp_e / (1-e*s(lat)**2)**0.5 # radius of curvature in the vertical prime
    alt = r*c(lat) + (z_p + e*N*s(lat))*s(lat) - N
    return [alt, lat, lon]

def latlongtoNED(H_LAN_LON):
    lon = H_LAN_LON[2]
    lat = H_LAN_LON[1]
    # Compute first in xyz coordinates(z: north pole, x - z plane: contains r, y: completes right - handed set)
    uDxyz = [-c(lat), 0, -s(lat)]
    uNxyz = [-s(lat), 0, c(lat)]
    uExyz = [0,1,0]

    # Rotate by longitude to change to PCPF frame
    L3 = [[c(lon), -s(lon), 0],
         [s(lon), c(lon), 0],
          [0, 0, 1]]
    uN = np.inner(L3, uNxyz)
    uE = np.inner(L3, uExyz)
    uD = np.inner(L3, uDxyz)

    return [uD, uN, uE]

# def rtolatlong(r_p,planet):
#     x_p = r_p.x
#     y_p = r_p.y
#     z_p = r_p.z
#
#     lam = atan2(y_p,x_p)
#
#     r = planet.Rp_p
#     R = planet.Rp_e
#
#     f = (R-r)/R # flattering
#
#     # Initialization
#     x_0 = (x_p**2 + y_p**2)**0.5
#     z_0 = z_p
#
#     A, B, C = r*z_0, 2*R*x_0, 2*(R**2 - r**2)
#
#
#     if z_0 != 0:
#         t0 = -(1 - f)*(x_0/z_0) + np.sign(z_0)*(1+ ((1-f)*(x_0/z_0))**2)**0.5
#         tkprev = t0
#         ind = 0
#         k = 0
#         while ind == 0:
#             k = k+1
#             f = A*tkprev**4 + (B+C)*tkprev**3 + (B-C)*tkprev -A
#             f_prime = 4*A*tkprev**3 + 3*(B+C)*tkprev**2 + (B-C)
#             tk = tkprev - f/f_prime
#             if np.abs(tk - tkprev) <= 10**-7:
#                 print('break')
#                 break
#             if k == 20:
#                 ind = 1
#             tkprev = tk
#         if ind == 0:
#             L = at(1/((1-f)*(1-tk**2)/(2*tk)))
#             h = np.sign(L)*(z_0 - r*(2*tk/(1+tk**2)))*(1+((1-f)*(1-tk**2)/2*tk)**2)**0.5
#         else:
#             pass
#             print('pass')
#
#     if z_0 == 0:
#         L = 0
#         h = x_0-R
#
#
#     return [H_LAN_LON(h, L, lam)]
#
#
# #
# def rtolatlong(r_p,planet):
#     # From PCPF (planet centered/planet fixed) to Geodetic
#     x_p = r_p.x
#     y_p = r_p.y
#     z_p = r_p.z
#
#     lam = atan2(y_p,x_p)
#
#     b = planet.Rp_p
#     a = planet.Rp_e
#     e = (1-(b/a)**2)**0.5
#     r = (x_p**2 + y_p**2)**0.5
#
#     E = (a*r + (a**2 - b**2))/(b*np.abs(z_p))
#     F = (a*r - (a**2 - b**2))/(b*np.abs(z_p))
#     P = (4/3)*(E*F + 1)
#     Q = 2*(E**2 -F**2)
#     D = P**3 + Q**2
#     vi = -((D)**0.5 + Q)**(1/3) + ((D)**0.5 - Q)**(1/3)
#     G = ((E**2 + vi)**0.5 + E)/2
#     t = (G**2 + (F - vi*G)/(2*G - E))**0.5 - G
#     phi = np.sign(z_p)*at(2*a*t/(b*(1-t**2)))
#     h = (r-a)*c(phi) + (np.abs(z_p) - b*t)* np.abs(s(phi))
#
#     return [H_LAN_LON(h, phi, lam)]




