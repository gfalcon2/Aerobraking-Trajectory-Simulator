#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
# python ABTS/setup.py build_ext --inplace
#!/bin/sh

import os
import sys
import io
import argparse , os
import os
import numpy as np
import random
import time as timecode
# import config

# insert at 1, 0 is the script path (or '' in REPL)
# export PYTHONPATH=`pwd`
sys.path.append(os.getcwd() + '/ABTS')
from simulation.run import run_analysis

parser = argparse.ArgumentParser()
# Miscellaneous Simulation
parser.add_argument('--results' , type=int , default=1 , help='# Generate csv file for results True=1, False=0')
parser.add_argument('--passresults' , type=int , default=1 , help='# Pass results as output True=1, False=0')
parser.add_argument('--print_res' , type=int , default=1 , help='# Print some lines True=1, False=0')
parser.add_argument('--directory_results' , type=str , default='' , help='# Directory where to save the results')
parser.add_argument('--directory_MarsGram' , type=str , default='' , help='# Directory where MarsGram is')
parser.add_argument('--MarsGram_version' , type=int , default=0 , help='# MarsGram x file to use', choices=range(0,150))
parser.add_argument('--montecarlo_analysis' , type=int , default=1 , help='# Generate csv file for Montecarlo results True=1, False=0')#parser.add_argument('--MarsGram_version_switch' , type=int , default=1 , help='# MarsGram x file to use', choices=[0,1,2,3])
parser.add_argument('--plot' , type=int , default=1 , help='# Generate png plots of results True=1, False=0')
parser.add_argument('--filename' , type=int , default=1 ,
                    help='# Filename with specifics of simulation, True =1, False=0')
parser.add_argument('--machine' , type=str , default='Cluster' ,
                    choices=['Laptop' , 'Cluster' , 'Aero' , 'Desktop_Home','Karnap_Laptop'])
parser.add_argument('--integrator' , type=str , default='Python' , choices=['Costumed' , 'Python'] ,
                    help='# Costumed customed integrator, Python python library integrator, only for drag passage, others phases use RK$%')

# Type of Missions
parser.add_argument('--type_of_mission' , type=str , default='Aerobraking Campaign' ,
                    choices=['Drag Passage' , 'Orbits' , 'Aerobraking Campaign'])
parser.add_argument('--number_of_orbits' , type=int , default=1000 , help='# number of aerobraking passage')

# Physical Models
parser.add_argument('--planet' , type=int , default=1 , help='# Earth = 0, Mars = 1, Venus = 2')
parser.add_argument('--planettime' , type=float , default=0 , help='Initial time of the mission, sec. Important for J2 effect and rotation of the planet')
parser.add_argument('--gravity_model' , type=str , default='Inverse Squared and J2 effect' ,
                    choices=['Constant' , 'Inverse Squared' , 'Inverse Squared and J2 effect'])
parser.add_argument('--density_model' , type=str , default='MarsGram' ,
                    choices=['Constant' , 'Exponential' , 'MarsGram'])
parser.add_argument('--wind' , type=int , default=1 ,
                    help='# wind calculation only if density model is MarsGram True=1, False=0')
parser.add_argument('--aerodynamic_model' , type=str , default='Mach-dependent' ,
                    choices=['Cd and Cl Constant' , 'Mach-dependent' , 'No-Ballistic flight with axial coefficient'],
                    help = '# "Mach-dependent" specific for spacecraft shape - "No-Ballistic flight" specific for blunted-cone shape')
parser.add_argument('--thermal_model' , type=str , default='Maxwellian Heat Transfer' ,
                    choices=['Maxwellian Heat Transfer' , 'Convective and Radiative'], help = '# "Maxwellian Heat Transfer" specific for spacecraft shape - "Convective and Radiative" specific for blunted-cone shape')

# Rates
parser.add_argument('--trajectory_rate' , type=float , default=100,
                    help='# rate at which the trajectory in drag passage integrate using RK4')
parser.add_argument('--flash1_rate' , type=float , default=3 , help='# rate at which Control Mode-1 is called')
parser.add_argument('--save_rate' , type=float , default=3 , help='# rate at which the data trajectory are saved')
parser.add_argument('--control_in_loop' , type=int , default=0 ,
                    help='# control in loop, control called during integration of trajectory, full state knowledge')
parser.add_argument('--flash2_through_integration' , type=int , default=0 ,
                    help='# integration of the equations of motion and lambda to define time switches and revaluation second time switch')

# Body
parser.add_argument('--body_shape' , type=str , default='Spacecraft' , choices=['Spacecraft' , 'Blunted Cone'])
parser.add_argument('--max_heat_rate' , type=float , default=0.1 ,
                    help='# Max heat rate the heat rate control will start to react to')
parser.add_argument('--max_heat_load' , type=float , default=30 ,
                    help='# Max heat load the heat load control will not be overcomed')
parser.add_argument('--dry_mass' , type=float , default=411.0 , help='# Initial dry mass of body in kg')
parser.add_argument('--prop_mass' , type=float , default=50.0 , help='# Initial propellant mass of body in kg')
parser.add_argument('--reflection_coefficient' , type=float , default=0.9 ,
                    help='# Diffuse reflection sigma =0, for specular reflection sigma = 1')
parser.add_argument('--thermal_accomodation_factor' , type=float , default=1 ,
                    help='# Thermal accomodation factor, Shaaf and Chambre')
parser.add_argument('--angle_of_attack' , type=float , default=90.0 , help='# max angle of attack of solar panels')

# Engine
parser.add_argument('--thrust' , type=float , default=4 , help='# Maximum magnitude thrust in N')

# Control Mode
parser.add_argument('--control_mode' , type=int , default=0,
                    help='# Use Rotative Solar Panels Control:  False=0, Only heat rate=1, Only heat load=2, Heat rate and Heat load = 3')  # change this
parser.add_argument('--security_mode' , type=int , default=1 ,
                    help='# Security mode that set the angle of attack to 0 deg if predicted heat load exceed heat load limit')
parser.add_argument('--second_switch_reevaluation' , type=int , default=1 ,
                    help='# Reevaluation of the second switch time when the time is closer to it')

# do not change
parser.add_argument('--heat_load_sol' , type=int , default=0, ###IMPORTANT
                    help='# heat load solution #leave it to 0 and change it only for control mode = 2:  Max energy depletaion=0, Min energy depletion=1, One switch max-min=2, One switch min-max = 3')  # change this

parser.add_argument('--thrust_control' , type=str , default='None' ,
                    choices=['None' , 'Aerobraking Maneuver' , 'Drag Passage Firing'])
parser.add_argument('--phi' , type=float , default=180 , help='# Thrust Angle, deg')
parser.add_argument('--delta_v' , type=float , default=0.1 , help='# Delta-v of Aerobraking Manuver,m/s')

parser.add_argument('--apoapsis_targeting' , type=int , default=0 ,
                    help='# Apoapsis Targeting Enabled')
parser.add_argument('--ra_fin_orbit' , type=float , default=25000*1e3 , help='# target final apoapsis for the orbit, km')

# Initial Conditions
parser.add_argument('--initial_condition_type' , type=int , default=0 ,
                    help='#  Initial Condition ra,hp = 0, Initial Condition v,gamma =1')

parser.add_argument('--ra_initial_a' , type=float , default=28523.95 * 1e3 ,
                    help='# Initial Apoapsis Radius for for-loop in m')
parser.add_argument('--ra_initial_b' , type=float , default=50000*1e3 ,
                    help='# Final Apoapsis Radius for for-loop in m')
parser.add_argument('--ra_step' , type=float , default=5000*1e33 , help='# Step Apoapsis Radius for for-loop in m')

parser.add_argument('--hp_initial_a' , type=float , default=94500 ,
                    help='# Initial Periapsis Altitude for for-loop in m')
parser.add_argument('--hp_initial_b' , type=float , default=159000 ,
                    help='# Final Periapsis Altitude for for-loop in m')
parser.add_argument('--hp_step' , type=float , default=1000000 , help='# Step Periapsis Radius for for-loop in m')

parser.add_argument('--v_initial_a' , type=float , default=3700 ,
                    help='# Initial Velocity (m/s) for for-loop if initial conditions are in v and gamma')
parser.add_argument('--v_initial_b' , type=float , default=5000 ,
                    help='# Final Velocity (m/s) for for-loop if initial conditions are in v and gamma')
parser.add_argument('--v_step' , type=float , default=100 ,
                    help='# Step Velocity (m/s) for for-loop if initial conditions are in v and gamma')

parser.add_argument('--gamma_initial_a' , type=float , default=2.5 ,
                    help='# Initial Gamma (deg) for for-loop if initial conditions are in v and gamma')
parser.add_argument('--gamma_initial_b' , type=float , default=7 ,
                    help='# Final Gamma (deg) for for-loop if initial conditions are in v and gamma')
parser.add_argument('--gamma_step' , type=float , default=0.5 ,
                    help='# Step Gamma (deg) for for-loop if initial conditions are in v and gamma')

parser.add_argument('--inclination' , type=float , default=93.6 , help='# Inclination Orbit, deg')  # 93.6
parser.add_argument('--omega' , type=float , default=0.0 , help='# AOP, deg')
parser.add_argument('--OMEGA' , type=float , default=0.0 , help='# RAAN, deg')

parser.add_argument('--EI' , type=float , default=160.0 , help='# Entry Interface, km')
parser.add_argument('--AE' , type=float , default=160.0 , help='# Atmospheric Exit, km')

parser.add_argument('--year' , type=int , default=2001 , help='# Mission year')
parser.add_argument('--month' , type=int , default=11 , help='# Mission month')
parser.add_argument('--day' , type=int , default=6 , help='# Mission day')
parser.add_argument('--hours' , type=int , default=8 , help='# Mission hours')
parser.add_argument('--minutes' , type=int , default=30 , help='# Mission minutes')
parser.add_argument('--secs' , type=int , default=0 , help='# Mission seconds')

# Final Conditions
parser.add_argument('--final_apoapsis' , type=float , default=4905.974818462152 * 10 ** 3 ,
                    help='# Final apoapsis radius if aerobraking campaign')

# MonteCarlo Simulations
parser.add_argument('--montecarlo' , type=int , default=0 , help='# Run Monte Carlo simulation True=1, False=0')
parser.add_argument('--initial_montecarlo_number' , type=int , default=1 , help='# Initial Monte Carlo sample number')
parser.add_argument('--montecarlo_size' , type=int , default=1000 , help='# number of Monte Carlo samples')

# MonteCarlo Perturbations
parser.add_argument('--CD_dispersion' , type=float , default=10 ,
                    help='# Max dispersion of CD for Uniform Distribution, %')
parser.add_argument('--CL_dispersion' , type=float , default=10 ,
                    help='# Max dispersion of CL for Uniform Distribution, %')
parser.add_argument('--rp_dispersion' , type=float , default=2.5 ,
                    help='# Max dispersion for initial vacuum periapsis radius following uniform distribution, km')
parser.add_argument('--ra_dispersion' , type=float , default=2.5 ,
                    help='# Max dispersion for initial apoapsis radius following uniform distribution, km')
parser.add_argument('--i_dispersion' , type=float , default=0.25 ,
                    help='# Max dispersion for initial inclination following uniform distribution, deg')
parser.add_argument('--raan_dispersion' , type=float , default=0.25 ,
                    help='# Max dispersion for initial right ascension of the ascending node following uniform distribution, deg')
parser.add_argument('--aop_dispersion' , type=float , default=0.25 ,
                    help='# Max dispersion for initial argument of periapsis following uniform distribution, deg')
parser.add_argument('--vi_dispersion' , type=float , default=0.025 ,
                    help='# Max dispersion for initial true anomaly following uniform distribution, deg')

# MonteCarlo Perturbation Guidance - Closed Form Solution (only for online)
parser.add_argument('--ra_dispersion_gnc' , type=float , default=0.25 ,
                    help='# Max dispersion for initial apoapsis radius used by gnc following uniform distribution, km')
parser.add_argument('--rp_dispersion_gnc' , type=float , default=0.25 ,
                    help='# Max dispersion for initial periapsis radius used by gnc following uniform distribution, km')
parser.add_argument('--i_dispersion_gnc' , type=float , default=0.025 ,
                    help='# Max dispersion for initial inclination used by gnc following uniform distribution, deg')
parser.add_argument('--raan_dispersion_gnc' , type=float , default=0.025 ,
                    help='# Max dispersion for initial right ascension of the ascending node used by gnc following uniform distribution, deg')
parser.add_argument('--aop_dispersion_gnc' , type=float , default=0.025 ,
                    help='# Max dispersion for initial argument of periapsis used by gnc following uniform distribution, deg')
parser.add_argument('--vi_dispersion_gnc' , type=float , default=0.0025 ,
                    help='# Max dispersion for initial true anomaly used by gnc following uniform distribution, deg')
# - online trajectory control (heat rate)
parser.add_argument('--rho_mudispersion_gnc' , type=float , default=0 ,
                    help='# Mean dispersion of rho for Gaussian Distribution, %')
parser.add_argument('--rho_sigmadispersion_gnc' , type=float , default=1 ,
                    help='# Std dispersion of rho for Gaussian Distribution, %')
parser.add_argument('--T_mudispersion_gnc' , type=float , default=0 ,
                    help='# Mean dispersion of T for Gaussian Distribution, %')
parser.add_argument('--T_sigmadispersion_gnc' , type=float , default=1 ,
                    help='# Std dispersion of T for Gaussian Distribution, %')
parser.add_argument('--S_mudispersion_gnc' , type=float , default=0 ,
                    help='# Mean dispersion of S for Gaussian Distribution, %')
parser.add_argument('--S_sigmadispersion_gnc' , type=float , default=1 ,
                    help='# Std dispersion of S for Gaussian Distribution, %')

parser.add_argument('--multiplicative_factor_heatload' , type=float , default=1 ,
                    help='# Multiplicative factor for heat rate prediction when calculated heat load')

parser.add_argument('--Odyssey_sim', type=int, default =0, help='# Simulate Odyssey Mission')

if __name__ == "__main__":
    args = parser.parse_known_args()[0]

    t = timecode.time()
    sol = run_analysis(args)
    if args.passresults:
        print('Ra initial = {0:.2f}km, Ra new = {1:.2f}km - Actual periapsis altitude = {2:.2f}km - Target Ra = {3:.2f}km - COMPUTATIONAL TIME {4:.2f}'.format(((sol.orientation.oe[0][0] * (1 + sol.orientation.oe[1][0])) * 1e-3),((sol.orientation.oe[0][-1] * (1 + sol.orientation.oe[1][-1])) * 1e-3),
           (min(sol.orientation.alt) * 1e-3), (args.final_apoapsis * 1e-3), (timecode.time() - t)))






