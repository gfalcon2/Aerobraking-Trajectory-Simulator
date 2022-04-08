#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

import subprocess
import time
import config
import os
import pandas as pd
import numpy as np

def MarsGram_online(model):
    # Change input file
    keyword_month =    '  MONTH    ='
    keyword_day =      '  MDAY     ='
    keyword_year =     '  MYEAR    ='
    keyword_hour =     '  IHR      ='
    keyword_min =      '  IMIN     ='
    keyword_sec =      '  SEC      ='
    keyword_points =   '  NPOS     ='
    keyword_inlat =    '  FLAT     ='
    keyword_inlon =    '  FLON     ='
    keyword_inhgt =    '  FHGT     ='
    keyword_deltime =  '  DELTIME  ='
    keyword_dellat =   '  DELLAT   ='
    keyword_dellon =   '  DELLON   ='
    keyword_delhgt =   '  DELHGT   ='
    keyword_numberMC = '  NMONTE   ='
    keyword_seeds =    '  NR1      ='

    version = model['Version']
    file_object = model['Directory']+'Code/inputstd/inputstd'+str(version)+'.txt'

    lines = []
    altitude = str(model['Initial Altitude'])
    latitude = str(model['Initial Latitude'])
    longitude = str(model['Initial Longitude'])

    timereal = model['Time Real']
    year = str(timereal.year)
    month = str(timereal.month)
    day = str(timereal.day)

    hour = str(timereal.hour)
    minute = str(timereal.minute)
    sec = str(timereal.second)

    point_number = str(model['Number of Points'])
    delalt = str(model['Delta Altitude'])
    dellon = str(model['Delta Longitude'])
    dellat = str(model['Delta Latitude'])
    deltime = str(model['Delta t'])
    mc_number = str(model['Monte Carlo'])
    mc_seeds = str(int(config.index_MonteCarlo))

    with open(file_object) as infile:
        for line in infile:
            if keyword_month in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + month)
                #print(line)
            elif keyword_day in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + day)
                #print(line)
            elif keyword_year in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + year)
                #print(line)
            elif keyword_hour in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + hour)
                #print(line)
            elif keyword_min in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + minute)
                #print(line)
            elif keyword_sec in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + sec)
                #print(line)
            elif keyword_points in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + str(point_number))
                #print(line)
            elif keyword_inlat in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + latitude)
                #print(line)
            elif keyword_inlon in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + longitude)
                #print(line)
            elif keyword_inhgt in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + altitude)
                #print(line)
            elif keyword_delhgt in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + delalt)
                #print(line)
            elif keyword_dellon in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + dellon)
                #print(line)
            elif keyword_dellat in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + dellat)
                #print(line)
            elif keyword_deltime in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + deltime)
                #print(line)
            elif keyword_numberMC in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + mc_number)
                #print(line)
            elif keyword_seeds in line:
                index = line.find('=')
                old = line[index + 1:-1]
                line = line.replace(str(old), ' ' + mc_seeds)
            lines.append(line)

    outfile = open(file_object, 'w')
    outfile.writelines(lines)

    infile.close()
    outfile.close()

    start_time = time.time()
    marsgram_app = model['Directory']+'Code/marsgram_M10_'+str(version)+'.x'
    output_dir = model['Directory']+'Code/output/OUTPUT'+str(version)+'.txt'
    args = (marsgram_app , '-d' ,output_dir)  # , "-c", "somefile.xml", "-d", "text.txt", "-r", "aString", "-f", "anotherString")


    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    output = popen.stdout.read()
    #print("--- MARSGram execution %s seconds ---" % (time.time() - start_time))

    # Read Output File
    file_object = model['Directory']+'Code/output/OUTPUT'+str(version)+'.txt'


    data = pd.read_csv(file_object,delim_whitespace = True,skiprows=[0], header=None, names=['Time','Alt','Lat','Lon','Dens','Temp','EW','EWTot','NW','NWTot','VWind','1','2','3','4','5','6','7','8','9','10','11','12','13','14','DensP'])
    data.head()#, engine='python'
    temp = [[], [], [], [], [], [], [], [], [], [], []]
    temp[0] = data['Alt'].tolist()
    temp[1] = data['Lat'].tolist()
    temp[2] = data['Lon'].tolist()
    temp[3] = data['Dens'].tolist()
    temp[4] = data['Temp'].tolist()
    temp[5] = data['EW'].tolist()
    temp[6] = data['EWTot'].tolist()
    temp[7] = data['NW'].tolist()
    temp[8] = data['NWTot'].tolist()
    temp[9] = data['VWind'].tolist()
    temp[10] = data['DensP'].tolist()


    keys = ['HgtMOLA','LatPC','LonW','Denkgm3','Temp','EWind','EWTot','NWind','NWTot','VWind','DensP']
    config.atmospheric_data = {keys[i]: temp[i] for i in range(len(keys))}

