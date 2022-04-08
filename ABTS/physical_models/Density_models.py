# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import math
import numpy as np
from physical_models.MarsGram import MarsGram_online
import config
from utils.Ref_system_conf import OE as orbitalelementconversion
from utils.Ref_system_conf import cartesian
import utils.Reference_system as Reference_system

import random
# Density
def density_constant ( h , p , OE=0 , lat=0 , lon=0 , timereal=0 , t0=0 , tf_prev = 0, montecarlo=0 , Wind=0, args= 0,version=[]):
    if config.drag_state == False:
        rho = 0.0
    else:
        rho = p.rho_ref
    T = temperature_linear(h , p)
    wind = [0 , 0 , 0]
    return rho , T , wind


def density_no ( h , p , OE=0 , lat=0 , lon=0 , timereal=0 , t0=0 , tf_prev = 0,montecarlo=0 , Wind=0, args= 0, version=[]):
    T = temperature_linear(h , p)
    wind = [0 , 0 , 0]
    return 0 , T , wind


def density_exp ( h , p , OE=0 , lat=0 , lon=0 , timereal=0 , t0=0 , tf_prev = 0, montecarlo=0 , Wind=0, args=0, version=[]):
    rho = p.rho_ref * np.exp(np.divide(p.h_ref - h , p.H))
    T = temperature_linear(h , p)
    wind = [0 , 0 , 0]
    return rho , T , wind


def marsgram(h, p, OE, lat, lon, timereal, t0, tf_prev, montecarlo, Wind, args, version=[]):
    if not version:
        version = args.MarsGram_version
    if config.drag_state == False:
        rho , T , wind = density_exp(h=h , p=p)
    else:
        rho_found = False
        first_step = True  # if the marsgram script is just been called

        # Preparation Initial Condition Data
        a , e , i , OMEGA , omega , vi = OE[0] , OE[1] , OE[2] , OE[3] , OE[4] , OE[5]
        lat , lon = math.degrees(lat) , math.degrees(lon)
        if lon < 0:
            lon = 360.0 + lon  # Assure always positive longitude values

        model = {'Monte Carlo': int(montecarlo) , 'Directory': args.directory_MarsGram , 'Version':version,'Wind': int(Wind) , 'Time Real': timereal ,
                 'Initial Longitude': lon , 'Initial Latitude': lat , 'Initial Altitude': h*1e-3}
        # print('MarsGram')
        while rho_found == False:

            if (not bool(config.atmospheric_data)) or (bool(
                    config.atmospheric_data) and first_step == False):  # if a table of densities is not already been created or if a table has been created but we are not at first time running the while loop
                # Create parameters:
                # descending in the atmosphere
                if (vi <= 2 * math.pi) and (vi > math.pi):
                    final_state_angle = 0.0
                    if not config.state_inner_boundary_atmosphere:  # save the initial right ascension for the first boundary of the atmosphere
                        config.state_inner_boundary_atmosphere = vi
                    model = define_modelparameters(OE , p , t0 , tf_prev, model , final_state_angle)
                else:
                    # ascending in the atmosphere
                    final_state_angle = np.pi *2 +0.05 - config.state_inner_boundary_atmosphere  # rad
                    model = define_modelparameters(OE , p , t0 , tf_prev, model , final_state_angle)

                if first_step == False:
                    config.previous_atmospheric_data = config.atmospheric_data

                # Create Table
                # print('Initial Altitude', model['Initial Altitude'], 'Final Altitude', model['Final Altitude'])
                MarsGram_online(model)
                config.MarsGram_justrecalled = 1
                config.index_Mars_Gram_call += 1

            ## Interpolate Results
            # get tables
            table = config.atmospheric_data
            Hgt = table.get('HgtMOLA')

            # look for boundaries
            error = [abs(float(x) - model['Initial Altitude']) for x in Hgt]
            try:
                if int(sorted(error)[0]*10) != 0:
                    # print('error1',sorted(error)[0],'Alt Bond',Hgt[0],Hgt[-1],'Alt',model['Initial Altitude'])
                    if bool(config.previous_atmospheric_data):
                        table = config.previous_atmospheric_data
                        Hgt = table.get('HgtMOLA')

                    # look for boundaries
                        error = [abs(float(x) - model['Initial Altitude']) for x in Hgt]
                        if int(sorted(error)[0] * 10) != 0:
                        # print('error2', sorted(error)[0], 'Alt Bond', Hgt[0], Hgt[-1], 'Alt',
                        #     model['Initial Altitude'])
                            first_step = False
                            continue
                    else:
                        first_step = False
                        continue
            except:
                continue

            if len(error) != 0:

                indexes = sorted(range(len(error)) , key=error.__getitem__)
                # boundary 1
                index_b1 = indexes[0]
                index_b2 = indexes[1]

                # interpolate
                if (Hgt[index_b2]) == (Hgt[index_b1]):
                    x = 0
                else:
                    x = (model['Initial Altitude'] - (Hgt[index_b1])) / ((Hgt[index_b2]) - (Hgt[index_b1]))

                # rhoa = table['Denkgm3'][index_b1]
                # rhob = table['Denkgm3'][index_b2]

                if model['Monte Carlo'] == 0:
                    rhoa = table['Denkgm3'][index_b1]
                    rhob = table['Denkgm3'][index_b2]
                    rho, T = interp(rhoa, rhob, x), interp(table['Temp'][index_b1], table['Temp'][index_b2], x)
                else:
                    rhoa = table['Denkgm3'][index_b1]*table['DensP'][index_b1]
                    rhob = table['Denkgm3'][index_b2]*table['DensP'][index_b2]
                    rho, T = interp(rhoa, rhob, x), interp(table['Temp'][index_b1], table['Temp'][index_b2], x)

                if model['Wind'] == 1:
                    wind = wind_def(model , indexes , x, table)
                else:
                    wind = [0.0 , 0.0 , 0.0]

                h_measured , latitude_measured , longitude_measured = interp((Hgt[index_b1]) , (Hgt[index_b2]) , x) , (
                    interp(table['LatPC'][index_b1] , table['LatPC'][index_b2] ,
                           x)) , (interp(table['LonW'][index_b1] ,
                                         table['LonW'][index_b2] , x))
                # ERROR! if interpolates between 360 and 0

                # Error in latitude and longitude:
                lat_error , lon_error , h_error = abs(lat - latitude_measured) , abs(lon - longitude_measured) , (
                    abs(model['Initial Altitude'] - h_measured))
                # to reset periodicity of the two angles
                if lon_error > 300.0:
                    lon_error = abs(lon_error - 360.0)

                first_step = False

                if ((((lat_error) + (lon_error)) < 100) and (abs(
                        x) <= 1.0)) or config.MarsGram_recall == 0:  # 7if latitude and longitude error less than a threshold and the point is inside the grid
                    rho_found = True
                    config.MarsGram_recall = 1
                    # if abs(x) > 1:
                    #    print('Help! Check the grid, point outside!')
                elif rho < 0.0:
                    rho_found = False

            else:
                print('HELP! SOMETHING BIG WRONG IN MARS GRAM!!')
    # print(rho,T,wind)
    if rho < 0.0:
        rho = 0.0
    return rho , T , wind


def temperature_linear ( h , p ):
    # into atmosphere
    if config.drag_state == True:
        T = p.T  # From Vikings data https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/JB095iB09p14811 (it should be higher)
    else:
        T = p.T  # p.T #K Outer space
    return T


def wind_def ( model , indexes , x, table):
    # boundary 1
    index_b1 = indexes[0]
    index_b2 = indexes[1]
    if model['Monte Carlo'] == 0:
        wE , wN = interp(table['EWind'][index_b1] ,
                         (table['EWind'][index_b2]) , x) , interp(
            table['NWind'][index_b1] , table['NWind'][index_b2] , x)

        wind = [wE , wN , 0]

    elif model['Monte Carlo'] == 1:
        wE , wN , vW = interp(table['EWTot'][index_b1] ,
                              (table['EWTot'][index_b2]) , x) , interp(
            table['NWTot'][index_b1] , table['NWTot'][index_b2] , x) , interp(
            table['VWind'][index_b1] , table['VWind'][index_b2] , x)
        wind = [wE , wN , vW]
    return wind


def define_modelparameters( OE , p , t0 , tf_prev, model , final_state_angle ):
    a , e , i , OMEGA , omega , vi = OE[0] , OE[1] , OE[2] , OE[3] , OE[4] , OE[5]
    # Create parameters:
    OMEGA_rate = (- (1.5 * ((p.mu) ** 0.5 * p.J2 * p.Rp_e ** 2) / ((1 - e ** 2) * a ** (7.0 / 2.0))) * math.cos(
        i))  # rad/s
    omega_rate = (OMEGA_rate * (((5.0 / 2.0) * (math.sin(i) ** 2) - 2.0) / math.cos(i)))  # rad/s

    # define initial and final state
    initial_state_angle = vi
    if e <1:
        E_initialstate = 2.0 * math.atan(
        ((1 - e) / (1 + e)) ** 0.5 * math.tan((initial_state_angle) *0.5))  # eccentric anomaly
        E_finalstate = 2.0 * math.atan(((1 - e) / (1 + e)) ** 0.5 * math.tan((final_state_angle)*0.5))  # eccentric anomaly

        # evaluate time to reach next state
        delta_t = abs((a ** 3.0 / p.mu) ** 0.5 * (
                (E_finalstate - e * math.sin(E_finalstate)) - (E_initialstate - e * math.sin(E_initialstate))))
    elif e>1:
        F_initialstate = 2.0*math.atanh(((e-1)/(e+1))**0.5*math.tan(vi/2))
        F_finalstate = 2.0*math.atanh(((e-1)/(e+1))**0.5*math.tan(0.0))
        delta_t =  abs((-a ** 3.0 / p.mu) ** 0.5 * ((e*math.sinh(F_finalstate)-F_finalstate)-(e*math.sinh(F_initialstate)-F_initialstate)))
        OMEGA_rate = 0
        omega_rate = 0
    # Update orbital parameters because of the Earth's oblateness
    omega_min = (omega + omega_rate * delta_t)
    OMEGA_min = (OMEGA + OMEGA_rate * delta_t)

    # Evaluate final state
    # print('final state angle', math.degrees(final_state_angle), 'vi', math.degrees(vi),'alt', model['Initial Altitude'])
    if final_state_angle > 0.0 and vi>final_state_angle:
        final_state_angle = vi+0.1
    # print('final state angle', math.degrees(final_state_angle), 'vi', math.degrees(vi),'alt', model['Initial Altitude'])
    oe = [a , e , i , OMEGA_min , omega_min , final_state_angle , 0]

    r , v = Reference_system.orbitalelemtorv(oe , p)
    r = np.array(r)

    # From PCI (planet centered inertial) to PCPF (planet centered/planet fixed)
    theta = (np.linalg.norm(p.omega) * (t0 + delta_t + tf_prev))

    R3 = np.array([[math.cos(theta) , math.sin(theta) , 0] , [-math.sin(theta) , math.cos(theta) , 0] , [0 , 0 , 1]])
    r_xprime = np.inner(R3 , r)

    LatLong = Reference_system.rtolatlong(r_xprime , p)

    model['Number of Points'] = int(abs(delta_t) * 3)
    if model['Number of Points'] == 0:
        model['Number of Points'] = 1

    model['Final Latitude'] = math.degrees(LatLong[1])
    model['Final Longitude'] = math.degrees(LatLong[2])
    model['Final Altitude'] = LatLong[0] * 1e-3

    if model['Final Longitude'] < 0:
        model['Final Longitude'] = 360 + model['Final Longitude']

    model['Delta Longitude'] = model['Final Longitude'] - model['Initial Longitude']

    if (model['Final Longitude'] - model['Initial Longitude']) > 250.0:  # if passage to 0 degree from final to initial
        model['Delta Longitude'] = 360.0 - model['Delta Longitude']
    elif (model['Final Longitude'] - model['Initial Longitude']) < -250.0:  # if passage to 0 degree from initial to final
        model['Delta Longitude'] = 360.0 + model['Delta Longitude']

    model['Delta Longitude'] = model['Delta Longitude'] / model['Number of Points']
    model['Delta Latitude'] = (model['Final Latitude'] - model['Initial Latitude']) / model['Number of Points']
    model['Delta Altitude'] = (model['Final Altitude'] - model['Initial Altitude']) / model[
        'Number of Points']
    model['Delta t'] = delta_t / model['Number of Points']

    model['Number of Points'] = model['Number of Points'] + 50
    return model


def interp ( a , b , x ):
    # check delta == diff b and a
    if (abs(b - a) > 20.0):
        if b <= (360.0) and b >= 350.0:
            b = 360.0 - b
        elif a <= (360.0) and a >= 350.0:
            a = 360.0 - a
    value = x * (b - a) + a
    return value
