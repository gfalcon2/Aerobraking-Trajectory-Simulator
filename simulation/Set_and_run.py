# Main file to define purposes, principal loop, input files, graph and results
from physical_models.Mission import *
from simulation.Aerobraking import aerobraking
import config as cnf
import math
from physical_models.Planet_data import planet_data
import numpy as np
import time as timecode
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import *
import os
from utils.save_cvs import *


def aeroraking_campaign(state,simulation):
    # Save Results??
    save_res = simulation['Results']


    # Descent towards Mars
    def mission_conf ():
        mission = {'Purpose': purpose , 'Gravity Model': simulation['Gravity Model'] , 'Density Model': simulation['Density Model'] , 'Wind': simulation['Wind'],
                   'Aerodynamic Model': simulation['Aerodynamic Model'] , 'Thermal Model': simulation['Thermal Model'], 'Control': simulation['Control Model'], 'Monte Carlo': simulation['Monte Carlo']}
        return (mission)


    # First phase
    purpose = 'Aerobraking around Mars.'


    mission = mission_conf()
    if simulation['Print'] == True:
        ('mission is', mission)

    ip = missiondef(mission_conf())
    p_class = planet_data(ip.M.planet)

    # Vehicle - calculation notebook page1
    dry_mass = 411 #kg
    prop_mass = 50 #kg
    mass = dry_mass + prop_mass
    area_body = 7.26#33.38#7.26# This is recalculated for the new sc config. 11 (look notes)# m^2 2001 Mars Odyssey Aerobraking, Smith & Bell paper
    length_sp = 3.7617#11.4#3.7617#5.7 # m solar array length https://www.jpl.nasa.gov/news/press_kits/odysseyarrival.pdf
    height_sp = area_body/length_sp
    # Main Body
    length_ody = 2.2 # m
    height_ody = 1.7 # m
    width_ody =  2.6 # m

    # Time
    year_i = state['Year']
    month_i = state['Month']
    day_i = state['Day']
    hour_i = state['Hour']
    min_i = state['Min']
    second_i = state['Second']



    apoapsis = state['Apoapsis']# 5000*10**3#7523.95142557378*10**3#28523.95142557378*10**3 # From phase 1 30371.25610549734*10**3#43090.01227466145 *10**3 #30371.25610549734*10**3#25000 *10**3 # from Phase 2 28523.95142557378*10**3 km
    periapsis = p_class.Rp_e + (state['Periapsis'])*10**3 #98
    #print('periapsis', periapsis/10**3-3396, 'km')

    if simulation['Monte Carlo'] == True:
        np.random.seed(int(cnf.index_MonteCarlo))
        r = np.random.uniform(-2.5*10**3, +2.5*10**3,2) # uncertanties +2.5 km
        s = np.random.uniform(-0.25, +0.25,3) # uncertanties +2.5 km
        apoapsis += r[0]
        periapsis += r[1]
        state['Inclination'] += s[0]
        state['OMEGA'] += s[1]
        state['omega'] += s[2]


    # Initial Condition Calcs
    semimajoraxis_in = (apoapsis + periapsis)/2
    eccentricity_in = (apoapsis - periapsis)/ (apoapsis + periapsis)
    state['vi'] = vi = np.radians(180.0001)

    # Initial Condition

    try:
        if 'Only DragPassage' in simulation.keys():
            if (simulation['Only DragPassage'] == True):
                #if (simulation['IC v-r'] == True):  # The initial conditions in orbital elements need to be calculated:
                r = p_class.Rp_e + 160 * 10 ** 3
                state['vi'] = - math.acos(
                        1 / eccentricity_in * (semimajoraxis_in * (1 - eccentricity_in ** 2) / r - 1))
    except:
        pass

    ## Initial Model Definition
    # Body
    def body ():
        Mass = mass  # kg
        length_SA = length_sp
        height_SA = height_sp
        Area_SA = length_SA*height_SA  # m^2
        length_SC = length_ody
        height_SC = height_ody
        Area_SC = length_ody*height_ody  # m^2
        b = cnf.model.body(Mass, length_SA, height_SA, Area_SA, length_SC, height_SC, Area_SC)
        return b
    b_class = body()


    # Initial Condition Entry
    def initialconditions ():
        a = semimajoraxis_in        #semi-major axis
        e = eccentricity_in         #eccentricity
        i = np.radians(state['Inclination'])                       #inclination # use 89.99 for 90
        OMEGA = np.radians(state['OMEGA'])                   #longitude of the ascending node
        omega = np.radians(state['omega'])                   #argument of the periapsis
        vi = state['vi']  # deg
        m0 = mass  # kg
        year = state['Year']
        month = state['Month']
        day = state['Day']
        hour = state['Hour']
        min = state['Min']
        second = state['Second']
        ic = cnf.model.initialcondition(a , e , i , OMEGA , omega, vi, m0, year, month, day, hour, min, second)
        return ic
    ic_class = initialconditions()

    # Aerodynamics
    def aerodynamics ():
        delta = math.radians(0)
        aoa = math.radians(90)
        thermal_accomodation_factor = 1.0 #0.2
        accomodation_factor = 0.9 # for diffuse reflection a_f =0, for specular reflection a_f = 1
        thermal_contact = 0 # for thermal perfect contact between rear and frontal area t_f = 0, for thermal insulation between rear and frontal area t_f = 1,
        thermal_limit = state['Heat Rate'] #W/cm^2
        a = cnf.model.aerodynamics(delta, aoa, thermal_accomodation_factor, accomodation_factor, thermal_contact, thermal_limit)
        return a
    a_class = aerodynamics()

    # Engine
    def engine ():
        phi = 0  # deg
        g_e = 9.81  # m/s
        T = 0  # N
        Isp = 300  # s
        e = cnf.model.engine(phi , g_e , T , Isp)
        return e
    e_class = engine()


    def model():
        body = b_class
        planet = p_class
        initialcondition = ic_class
        aerodynamics = a_class
        engine = e_class
        m = cnf.model(body, planet, initialcondition, aerodynamics, engine)
        return m
    m = model()

    # Define Final Condition
    Period_final = 2 #h
    index,semimajoraxis_fin = 0, state['Final SMA'] #End of third phase (2000+3396)*10**3#(((Period_final*3600)/(2*math.pi))**2 * m.planet.mu)**(1/3)
    terminal_state = [index, semimajoraxis_fin]

    # RUN SIMULATION
    t = timecode.time()
    aerobraking(ip, m, terminal_state, simulation)
    elapsed = timecode.time() - t

    print('rho: ',max(cnf.solution.physical_properties.rho), 'kg/m^3')
    print('heat rate: ',max(cnf.solution.performance.heat_rate), 'W/cm^2')

    # Save results
    if save_res == 1:
        now = datetime.now()
        if 'Filename' in simulation.keys():
            if simulation['Monte Carlo'] == True:
                folder_name = simulation['Filename'][0:simulation['Filename'].find('_nMC')]
            else:
                folder_name= simulation['Filename']
            name = '/Users/giusyfalcone/Aerobraking_SA_project_results/' + folder_name + '/' + simulation['Filename']
            filename = name + '.csv'
        else:
            name = '/Users/giusyfalcone/Aerobraking_SA_project_results/' + str(date.today()) + '/' + str(now.strftime("%H_%M_%S"))
            filename = name + '.csv'
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        save_cvs(filename)

    print('Elapsed Time :', elapsed)

    if simulation['Plot'] == True:
        ## Plot

        fig = plt.figure()
        fig.set_size_inches(20.5, 19)
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(cnf.solution.orientation.pos_ii[0],cnf.solution.orientation.pos_ii[1], cnf.solution.orientation.pos_ii[2])
        wframe = None
        tstart = timecode.time()

        oldcol = wframe
        r = random.randint(1, m.planet.Rp_e)
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = r * np.outer(np.cos(u), np.sin(v))
        y = r * np.outer(np.sin(u), np.sin(v))
        z = r * np.outer(np.ones(np.size(u)), np.cos(v))
        sphere = ax.plot_surface(x, y, z, color='b')

        ax.set_xlabel('X')
        ax.set_xlim(-apoapsis, periapsis)
        ax.set_ylabel('Y')
        ax.set_ylim(-apoapsis, periapsis)
        ax.set_zlabel('Z')
        ax.set_zlim(-m.planet.Rp_e, m.planet.Rp_e)
        ax.view_init(220, 255)
        fig.savefig(name + 'figure1.png')



        circle1 = plt.Circle((0, 0), m.planet.Rp_e, color='r')
        fig2 = plt.figure()

        #ax = fig.add_subplot(111, projection='3d')
        plt.plot(cnf.solution.orientation.pos_ii[0],cnf.solution.orientation.pos_ii[1])
        plt.plot(cnf.solution.orientation.pos_ii[0][0],cnf.solution.orientation.pos_ii[1][0],'.',c='red')
        #fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
        # (or if you have an existing figure)
        fig = plt.gcf()
        fig.set_size_inches(20.5, 4)
        ax = fig.gca()
        ax.add_artist(circle1)
        fig2.savefig(name + 'figure2.png')

        fig3 = plt.figure()
        plt.plot([item/(60*60) for item in cnf.solution.orientation.time], cnf.solution.forces.energy,'.')
        plt.xlabel('Time, hr')
        plt.ylabel('Energy, J/kg')
        fig3.savefig(name + 'figure3.png')



        fig4 = plt.figure()
        plt.plot([item/(60*60) for item in cnf.solution.orientation.time], cnf.solution.performance.heat_rate)
        plt.xlabel('Time, hr')
        plt.ylabel('Heat Rate, W/cm^2')
        fig4.savefig(name + 'figure4.png')

        fig5 = plt.figure()
        plt.plot([item/(60*60) for item in cnf.solution.orientation.time], cnf.solution.performance.heat_load)
        plt.xlabel('Time, hr')
        plt.ylabel('Heat, J/cm^2')
        fig5.savefig(name + 'figure5.png')

        fig6 = plt.figure()
        plt.plot([item/(10**3) for item in cnf.solution.orientation.alt], cnf.solution.performance.heat_rate)
        plt.xlabel('Altitude, km')
        plt.ylabel('Heat Rate, W/cm^2')
        fig6.savefig(name + 'figure6.png')

        fig7 = plt.figure()
        plt.plot([item/(60*60) for item in cnf.solution.orientation.time], [math.degrees(item2) for item2 in cnf.solution.physical_properties.aoa])
        plt.xlabel('Time, s')
        plt.ylabel('Angle of Attack, degrees')
        fig7.savefig(name + 'figure7.png')
        plt.close()
