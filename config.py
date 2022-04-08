#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""

import numpy as np

class model:
    def __init__(self, body, planet,initialcondition, aerodynamics, engine):
        self.body = body
        self.planet = planet
        self.initialcondition = initialcondition
        self.aerodynamics = aerodynamics
        self.engine = engine

    # Body
    class body:
        def __init__(self, Mass=0, length_SA=0, height_SA=0, Area_SA=0, length_SC=0, height_SC=0, Area_SC=0, Area_tot = 0, delta=0, NoseRadius=0, BaseRadius=0):
            self.Mass = Mass
            self.length_SA = length_SA
            self.height_SA = height_SA
            self.Area_SA = Area_SA
            self.length_SC = length_SC
            self.height_SC = height_SC
            self.Area_SC = Area_SC
            self.Area_tot = Area_tot
            self.delta = delta
            self.NoseRadius = NoseRadius
            self.BaseRadius = BaseRadius

#    def MCMass(self, Mass):
#        self.Mass = Mass
    # Planet
    class planet:
        def __init__(self, Rp_e, Rp_p, Rp_m, mass, p, k, omega, g_ref, rho_ref, h_ref, H, R, gamma, T, J2, mu, mu_fluid, Lz):
            self.Rp_e = Rp_e
            self.Rp_p = Rp_p
            self.Rp_m = Rp_m
            self.mass = mass
            self.p = p
            self.k = k
            self.omega = omega
            self.g_ref = g_ref
            self.rho_ref = rho_ref
            self.h_ref = h_ref
            self.H = H
            self.R = R
            self.gamma = gamma
            self.T = T
            self.J2 = J2
            self.mu = mu
            self.mu_fluid = mu_fluid
            self.Lz = Lz


    # Aerodynamics
    class aerodynamics:
        def __init__(self, delta, aoa, thermal_accomodation_factor, reflection_coefficient, thermal_contact, heat_rate_limit,heat_load_limit):
            self.delta = delta
            self.aoa = aoa
            self.thermal_accomodation_factor = thermal_accomodation_factor
            self.reflection_coefficient = reflection_coefficient
            self.thermal_contact = thermal_contact
            self.heat_rate_limit = heat_rate_limit
            self.heat_load_limit = heat_load_limit


    class initialcondition:
        def __init__(self, a, e, i, OMEGA, omega, vi, m, year, month, day, hour, min, second, time_rot):
            self.a = a
            self.e = e
            self.i = i
            self.OMEGA = OMEGA
            self.omega = omega
            self.vi = vi
            self.m = m
            self.year = year
            self.month = month
            self.day = day
            self.hour = hour
            self.min = min
            self.second = second
            self.time_rot = time_rot

    # Engine
    class engine:
        def __init__(self, phi, g_e, T, Isp):
            self.phi = phi
            self.g_e = g_e
            self.T = T
            self.Isp = Isp


# Solution
impact = False
altitudeperiapsis = []
max_heatrate = []
solution_intermediate = []#np.empty((77,1))
atmospheric_data = {}
previous_atmospheric_data = {}
drag_state = False
ascending_phase = False
evaluate_switch_heat_load = False
security_mode = False
time_IEI = 0
time_OEI = 0
time_switch_1 = 0 # initialization of the first switch time for control 2 and 3
time_switch_2 = 0 # initialization of the second switch time for control 2 and 3
state_inner_boundary_atmosphere = []
count_aerobraking = 0 # Counter all aerobraking
count_dori = 0 # Counter for one passage
count_phase = 0
count_numberofpassage = 0
count_overcome_hr = 0
counter_random = 0
save_index_heat = 0
index_warning_alt = 0
index_warning_flow = 0
index_Mars_Gram_call = 0
index_MonteCarlo = 1#1024
index_propellant_mass = 1
T_w= 4
delta_v_man = 0
closed_form_solution_off = 1 # closed form solution never been called online by the GNC for a new passage
aoa = np.pi/2
aoa_past = np.pi/2
#Results to delete
periapsis_list= []
delta_v_list = []
orbit_number_list = []
heat_load_past = 0 #list used for flash1and3. List collecting all heat rate already experienced
heat_load_ppast = 0 #list used for flash1and3. Two steps before
state_flesh1 = []
aoa_list = []
initial_position_closed_form = []
continue_simulation = True # to continue to run second step integrator
timer_revaluation = 0
MarsGram_recall = 0
heat_rate_prev = 0
sensible_loads = False
counter_integrator =0
prev_step_integrator = 0
initial_time_save = 0

class ctrller():
    def __init__(self):
        self.guidance_t_eval = []
        self.count_controller = 1
        self.count_prev_controller = 0
        self.stored_state = 1
        self.prev_time = 0
        self.t = 0

class solution:
    def __init__(self, orientation, physical_prop, performance, forces, initial_condition,simulation):
        self.orientation = orientation
        self.physical_prop = physical_prop
        self.performance = performance
        self.forces = forces
        self.initial_condition = initial_condition
        self.simulation = simulation

    class orientation:
        def __init__(self, *args, **kwargs):
            self.time = []
            self.year = []
            self.month = []
            self.day = []
            self.hour = []
            self.min = []
            self.second = []
            self.numberofpassage = []
            #self.time.append(t)
            self.pos_ii = [[],[],[]]
            self.vel_ii = [[],[],[]]
            self.pos_ii_mag = []
            self.vel_ii_mag = []
            self.pos_pp = [[],[],[]]
            self.pos_pp_mag = []
            self.vel_pp = [[],[],[]]
            self.vel_pp_mag = []
            self.oe = [[], [], [], [], [], []]
            self.lat = []
            self.lon = []
            self.alt = []
            self.gamma_ii = []
            self.gamma_pp = []
            self.h_ii = [[], [], []]
            self.h_pp = [[],[],[]]
            self.h_ii_mag = []
            self.h_pp_mag = []
            self.uD = [[],[],[]]
            self.uE = [[],[],[]]
            self.uN = [[],[],[]]
            self.vN = []
            self.vE = []
            self.azi_pp = []

        def add(self, x):
            self.time.append(x)


    class physical_properties:
        def __init__(self):
            self.rho = []
            self.T = []
            self.p = []
            self.wind = [[], [], []]
            self.cL = []
            self.cD = []
            self.aoa = []
            self.S = []

    class performance:
        def __init__(self):
            self.mass = []
            self.heat_rate = []
            self.heat_load = []
            self.T_r = []
            self.q = []

    class forces:
        def __init__(self):
            self.gravity_ii = [[], [], []]
            self.drag_pp = [[], [], []]
            self.drag_ii = [[], [], []]
            self.lift_pp = [[], [], []]
            self.lift_ii = [[], [], []]
            self.force_ii = [[], [], []]
            self.energy = []
    class simulation:
        def __init__(self):
            self.MC_seed = []
            self.drag_passage = []
    class closed_form:
        def __init__(self):
            self.t_cf = []
            self.h_cf = []
            self.gamma_cf = []
            self.v_cf = []
solution.orientation = solution.orientation()
solution.physical_properties = solution.physical_properties()
solution.performance = solution.performance()
solution.forces = solution.forces()
solution.simulation = solution.simulation()
solution.closed_form = solution.closed_form()


# Performance
#    class performance:
#        def __init__(self, acceleration_peak, heatrate_peak, heatload_peak, deployaltitude):
#            self.acceleration_peak = acceleration_peak
#            self.heatrate_peak = heatrate_peak
#            self.heatload_peak = heatload_peak
#            self.deployaltitude = deployaltitude
