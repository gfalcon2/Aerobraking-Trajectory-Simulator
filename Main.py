#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

import os
import sys
import io
import argparse, os
from simulation.Run_analysis import run_analysis

parser = argparse.ArgumentParser()
# Miscellaneous Simulation
parser.add_argument('--results', type=int, default=1, help='# Generate csv file for results True=1, False=0')
parser.add_argument('--print_res', type=int, default=1, help='# Print some lines True=1, False=0')
parser.add_argument('--directory_results',type=str,default='', help='# Directory where to save the results')
parser.add_argument('--directory_MarsGram',type=str,default='', help='# Directory where MarsGram is')
parser.add_argument('--plot', type=int, default=0, help='# Generate png plots of results True=1, False=0')
parser.add_argument('--filename', type=int, default=1, help='# Filename with specifics of simulation, True =1, False=0')

# Physical Models
parser.add_argument('--planet', type=int, default=1, help='# Earth = 0, Mars = 1, Venus = 2')
parser.add_argument('--gravity_model', type=str, default='Inverse Squared and J2 effect', choices=['Constant', 'Inverse Squared', 'Inverse Squared and J2 effect'])
parser.add_argument('--density_model', type=str, default='MARSGram', choices=['Constant', 'Exponential', 'MARSGram'])
parser.add_argument('--wind', type=int, default=1, help='# wind calculation only if density model is MARSGram True=1, False=0')
parser.add_argument('--aerodynamic_model', type=str, default='Mach-dependent', choices=['Cd and Cl Constant', 'Mach-dependent','No-Ballistic flight with axial coefficient'])
parser.add_argument('--thermal_model', type=str, default='Maxwellian Heat Transfer', choices = ['Maxwellian Heat Transfer','Convective','Radiative'])

# Body                   
parser.add_argument('--body_shape',type=str, default='Spacecraft', choices=['Spacecraft','Blunted Cone'])
parser.add_argument('--max_heat_rate', type=float, default=0.1,help='# Max heat rate the heat rate control will start to react to')
parser.add_argument('--max_heat_load', type=float, default=30,help='# Max heat load the heat load control will not be overcomed')
parser.add_argument('--dry_mass',type=float, default=411,help='# Initial dry mass of body in kg')
parser.add_argument('--prop_mass',type=float, default=50,help='# Initial propellant mass of body in kg')
parser.add_argument('--accomodation_factor', type=float, default=0.9, help='# Diffuse reflection a_f =0, for specular reflection a_f = 1')
parser.add_argument('--thermal_accomodation_factor', type=float, default=1, help='# Thermal accomodation factor, Shaaf and Chambre')
parser.add_argument('--angle_of_attack', type=float, default=45, help='# initial angle of attack of the body')

# Engine
parser.add_argument('--thrust',type=float, default=20,help='# Maximum magnitude thrust in N')
parser.add_argument('--phi',type=float, default=0,help='# Thrust Angle, deg')

# Control Mode
parser.add_argument('--control_mode', type=int, default=0, help='# Use Rotative Solar Panels Control:  False=0, Only heat rate=1, Only heat load=2, Heat rate and Heat load = 3') #change this
parser.add_argument('--thrust_control', type=str, default='None', choices = ['None','Aerobraking Maneuver','Drag Passage Firing'] )
parser.add_argument('--firing_time', type=float, default=10, help = '# Duration of Aerobraking Manuver,s' )

# Type of Missions
parser.add_argument('--type_of_mission', type=str, default = 'Drag Passage', choices=['Drag Passage','Orbits','Aerobraking Campaign'])
parser.add_argument('--number_of_orbits', type=int, default=1000, help='# number of aerobraking passage')
  
# Initial Conditions
parser.add_argument('--initial_condition_type', type=int, default=0, help='#  Initial Condition ra,hp = 0, Initial Condition v,gamma =1')

parser.add_argument('--ra_initial_a', type = float, default=28523.95*10**3, help='# Initial Apoapsis Radius for for-loop in m')
parser.add_argument('--ra_initial_b', type = float, default=30000*10**3, help='# Final Apoapsis Radius for for-loop in m')
parser.add_argument('--ra_step', type = float, default=5000*10**3, help='# Step Apoapsis Radius for for-loop in m')

parser.add_argument('--hp_initial_a', type = float, default=94500, help='# Initial Periapsis Altitude for for-loop in m')
parser.add_argument('--hp_initial_b', type = float, default=120000, help='# Final Periapsis Altitude for for-loop in m')
parser.add_argument('--hp_step', type = float, default=50000, help='# Step Periapsis Radius for for-loop in m')

parser.add_argument('--v_initial_a', type = float, default=3700, help='# Initial Velocity (m/s) for for-loop if initial conditions are in v and gamma')
parser.add_argument('--v_initial_b', type = float, default=5000, help='# Final Velocity (m/s) for for-loop if initial conditions are in v and gamma')
parser.add_argument('--v_step', type = float, default=100, help='# Step Velocity (m/s) for for-loop if initial conditions are in v and gamma')

parser.add_argument('--gamma_initial_a', type = float, default=2.5, help='# Initial Gamma (deg) for for-loop if initial conditions are in v and gamma')
parser.add_argument('--gamma_initial_b', type = float, default=7, help='# Final Gamma (deg) for for-loop if initial conditions are in v and gamma')
parser.add_argument('--gamma_step', type = float, default=0.5, help='# Step Gamma (deg) for for-loop if initial conditions are in v and gamma')
                    
parser.add_argument('--inclination', type = float, default=93.6, help='# Inclination Orbit')
parser.add_argument('--omega', type = float, default=0, help='# AOP and RAAN')

parser.add_argument('--year', type = int, default =2001, help='# Mission year')
parser.add_argument('--month', type = int, default =11, help='# Mission month')
parser.add_argument('--day', type = int, default =6, help='# Mission day')
parser.add_argument('--hours', type = int, default =8, help='# Mission hours')
parser.add_argument('--minutes', type = int, default =30, help='# Mission minutes')
parser.add_argument('--secs', type = int, default =0, help='# Mission seconds')

# Final Conditions
parser.add_argument('--final_apoapsis', type=float,default=4905.974818462152 * 10 ** 3, help='# Final apoapsis radius if aerobraking campaign')

# MonteCarlo Simulations
parser.add_argument('--montecarlo', type=int, default=0, help='# Run Monte Carlo simulation True=1, False=0')
parser.add_argument('--montecarlo_size', type=int, default=100, help='# number of Monte Carlo samples')

# MonteCarlo Perturbations
parser.add_argument('--CD_dispersion',type=float,default=10,help='# Max dispersion of CD for Uniform Distribution, %')
parser.add_argument('--CL_dispersion',type=float,default=10,help='# Max dispersion of CL for Uniform Distribution, %')
parser.add_argument('--rp_dispersion',type=float,default=2.5,help='# Max dispersion for initial vacuum periapsis radius following uniform distribution, km')
parser.add_argument('--ra_dispersion',type=float,default=2.5,help='# Max dispersion for initial apoapsis radius following uniform distribution, km')
parser.add_argument('--i_dispersion',type=float,default=0.25,help='# Max dispersion for initial inclination following uniform distribution, deg')
parser.add_argument('--raan_dispersion',type=float,default=0.25,help='# Max dispersion for initial right ascension of the ascending node following uniform distribution, deg')
parser.add_argument('--aop_dispersion',type=float,default=0.25,help='# Max dispersion for initial argument of periapsis following uniform distribution, deg')

if __name__ == "__main__":
    args = parser.parse_known_args()[0]
    run_analysis(args)





