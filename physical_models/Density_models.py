# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import math
import numpy as np
from physical_models.MARSGram import MarsGram_online
import config
from utils.Ref_system_conf import OE as orbitalelementconversion
from utils.Ref_system_conf import cartesian
import utils.Reference_system as Reference_system
# Density
def density_constant (h, p, OE=0, lat=0, lon=0, timereal=0, t0=0, montecarlo=0, directory = 0,Wind=0) :
    rho = p.rho_ref
    T = temperature_linear(h, p)
    wind = [0 , 0 , 0]
    return rho, T, wind

def density_no ( h, p , OE=0, lat=0, lon=0, timereal=0,t0=0, montecarlo=0, directory=0,Wind=0) :
    T = temperature_linear(h, p)
    wind = [0 , 0 , 0]
    return 0, T, wind

def density_exp(h, p, OE=0, lat=0, lon=0, timereal=0,t0=0, montecarlo=0, directory = 0, Wind=0):
    rho = p.rho_ref * np.exp(np.divide(p.h_ref - h, p.H))
    T = temperature_linear(h, p)
    wind = [0 , 0 , 0]
    return rho, T, wind

def marsgram(h, p, OE, lat, lon, timereal ,t0, montecarlo, directory,Wind):
    if config.drag_state == False:
        rho,T,wind = density_exp(h=h, p=p)
    else:
        rho_found = False
        first_step = True  # if the marsgram script is just been called
        
        # Preparation Initial Condition Data
        lat, lon = math.degrees(lat), math.degrees(lon)
        if lon < 0:
            lon = 360 + lon # Assure always positive longitude values

        model = {'Monte Carlo': int(montecarlo),'Directory':directory, 'Wind': int(Wind), 'Time Real':timereal, 'Initial Longitude':lon, 'Initial Latitude': lat, 'Initial Altitude': h * 10 ** -3 }
        while rho_found == False:

            if (not bool(config.atmospheric_data)) or (bool(config.atmospheric_data) and first_step == False):  #if a table of densities is not already been created or if a table has been created but we are not at first time running the while loop
                # Create parameters:
                # descending in the atmosphere
                if (OE.vi<= 2*math.pi) and (OE.vi> math.pi):
                    final_state_angle = 0
                    if not config.state_inner_boundary_atmosphere: # save the initial right ascension for the first boundary of the atmosphere
                        config.state_inner_boundary_atmosphere = OE.vi
                    model = define_modelparameters(OE, p, t0, model, final_state_angle)
                else:
                    # ascending in the atmosphere
                    final_state_angle =6.4 - config.state_inner_boundary_atmosphere # rad
                    model = define_modelparameters(OE, p, t0, model, final_state_angle)


                # Create Table
                MarsGram_online(model)
                config.index_Mars_Gram_call += 1

            ## Interpolate Results
            # get tables
            Hgt = (config.atmospheric_data.get('HgtSFCM') or config.atmospheric_data.get('HgtMOLA'))

            # look for boundaries
            error = [abs(float(x) - model['Initial Altitude']) for x in Hgt]
            if len(error) != 0:

                indexes = sorted(range(len(error)), key=error.__getitem__)
                # boundary 1
                index_b1 = indexes[0]
                index_b2 = indexes[1]

                # interpolate
                if (Hgt[index_b2]) == (Hgt[index_b1]):
                    x = 0
                else:
                    x = (model['Initial Altitude'] - (Hgt[index_b1])) / ((Hgt[index_b2]) - (Hgt[index_b1]))

                rho,T = interp(config.atmospheric_data['Denkgm3'][index_b1], (config.atmospheric_data['Denkgm3'][index_b2]), x),  interp(config.atmospheric_data['Temp'][index_b1], config.atmospheric_data['Temp'][index_b2], x)

                if model['Wind'] == 1:
                    wind = wind_def(model, indexes, x)
                else:
                    wind = [0,0,0]

                h_measured,latitude_measured,longitude_measured = interp((Hgt[index_b1]),(Hgt[index_b2]),x),(interp(config.atmospheric_data['LatPC'][index_b1], config.atmospheric_data['LatPC'][index_b2], x)),  (interp(config.atmospheric_data['LonW'][index_b1], config.atmospheric_data['LonW'][index_b2], x))
                # ERROR! if interpolates between 360 and 0

                # Error in latitude and longitude:
                lat_error, lon_error, h_error = abs(lat-latitude_measured), abs(lon-longitude_measured), (abs(model['Initial Altitude'] - h_measured))
                # to reset periodicity of the two angles
                if lon_error > 300:
                    lon_error = abs(lon_error - 360)


                first_step = False

                if (((lat_error) + (lon_error)) < 10) and (abs(x)<=1): # if latitude and longitude error less than a threshold and the point is inside the grid
                    rho_found = True
                    if abs(x) > 1:
                        print('Help! Check the grid, point outside!')
            else:
                print('HELP! SOMETHING BIG WRONG IN MARS GRAM!!')

    return rho, T, wind


def temperature_linear(h, p):
    #into atmosphere
    if config.drag_state==True:
        T = p.T # From Vikings data https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/JB095iB09p14811 (it should be higher)
    else:
        T = p.T#p.T #K Outer space
    return T

def wind_def(model, indexes, x):
    # boundary 1
    index_b1 = indexes[0]
    index_b2 = indexes[1]
    if model['Monte Carlo'] == 0:
        wE, wN = interp(config.atmospheric_data['EWind'][index_b1],
                        (config.atmospheric_data['EWind'][index_b2]), x), interp(
            config.atmospheric_data['NWind'][index_b1], config.atmospheric_data['NWind'][index_b2], x)

        wind = [wE, wN, 0]

    elif model['Monte Carlo'] == 1:
        wE, wN, vW = interp(config.atmospheric_data['EWTot'][index_b1],
                            (config.atmospheric_data['EWTot'][index_b2]), x), interp(
            config.atmospheric_data['NWTot'][index_b1], config.atmospheric_data['NWTot'][index_b2], x), interp(
            config.atmospheric_data['VWind'][index_b1], config.atmospheric_data['VWind'][index_b2], x)
        wind = [wE, wN, vW]
    return wind


def define_modelparameters(OE, p, t0, model, final_state_angle):
    # Create parameters:
    OMEGA_rate = (- (3 / 2 * ((p.mu) ** 0.5 * p.J2 * p.Rp_e ** 2) / ((1 - OE.e ** 2) * OE.a ** (7 / 2))) * math.cos(
        OE.i))  # rad/s
    omega_rate = (OMEGA_rate * ((5 / 2 * (math.sin(OE.i) ** 2) - 2) / math.cos(OE.i)))  # rad/s


    # define initial and final state
    initial_state_angle = OE.vi
    E_initialstate = 2 * math.atan(((1 - OE.e) / (1 + OE.e)) ** 0.5 * math.tan((initial_state_angle) / 2)) # eccentric anomaly
    E_finalstate  = 2 * math.atan(((1 - OE.e) / (1 + OE.e)) ** 0.5 * math.tan((final_state_angle) / 2)) # eccentric anomaly

    # evaluate time to reach next state
    delta_t = (OE.a ** 3 / p.mu) ** 0.5 * ((E_finalstate - OE.e * math.sin(E_finalstate)) - (E_initialstate - OE.e * math.sin(E_initialstate)))

    # Update orbital parameters because of the Earth's oblateness
    omega_min = (OE.omega + omega_rate * delta_t)
    OMEGA_min = (OE.OMEGA + OMEGA_rate * delta_t)

    # Evaluate final state
    oe = orbitalelementconversion(OE.a, OE.e, OE.i, OMEGA_min, omega_min, final_state_angle, mass=0)
    r, v = Reference_system.orbitalelemtorv(oe, p)
    r = np.array([r.x, r.y, r.z])

    # From PCI (planet centered inertial) to PCPF (planet centered/planet fixed)
    theta = (np.linalg.norm(p.omega) * (t0 + delta_t))

    R3 = np.array([[math.cos(theta), math.sin(theta), 0], [-math.sin(theta), math.cos(theta), 0], [0, 0, 1]])
    r_xprime = np.inner(R3, r)
    r_xprime = cartesian(r_xprime[0], r_xprime[1], r_xprime[2])

    [LatLong] = Reference_system.rtolatlong(r_xprime, p)

    model['Number of Points'] = int(abs(delta_t) * 3)

    model['Final Latitude'] = math.degrees(LatLong.LAT)
    model['Final Longitude'] = math.degrees(LatLong.LON)
    model['Final Altitude'] = LatLong.h * 10 ** -3

    if model['Final Longitude'] < 0:
        model['Final Longitude'] = 360 + model['Final Longitude']

    model['Delta Longitude'] = model['Final Longitude'] - model['Initial Longitude']

    if (model['Final Longitude'] - model['Initial Longitude']) > 250:  # if passage to 0 degree from final to initial
        model['Delta Longitude'] = 360 - model['Delta Longitude']
    elif (model['Final Longitude'] - model['Initial Longitude']) < -250:  # if passage to 0 degree from initial to final
        model['Delta Longitude'] = 360 + model['Delta Longitude']


    model['Delta Longitude'] = model['Delta Longitude']/model['Number of Points']
    model['Delta Latitude'] = (model['Final Latitude'] - model['Initial Latitude']) / model['Number of Points']
    model['Delta Altitude'] = (model['Final Altitude'] - model['Initial Altitude'])  / model[
        'Number of Points']
    model['Delta t'] = delta_t / model['Number of Points']

    model['Number of Points'] = model['Number of Points'] + 50
    return model

def interp(a,b,x):
    # check delta == diff b and a
    if (abs(b-a) > 20):
        if b<= (360) and b >= 350:
            b = 360 - b
        elif a<= (360) and a >= 350:
            a = 360 - a
    value = x * (b - a) + a
    return value

