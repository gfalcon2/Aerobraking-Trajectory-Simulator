#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:23:54 2020

@author: Giusy Falcone (gfalcon2@illinois.edu)
@copyright University of illinois at Urbana Champaign
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec
import config as cnf
import numpy as np
from utils.misc import find_neighbour


def plots(state,m, name,args):
    traj_3d(state,m,name)
    traj_2d(state,m, name)
    #if not (cnf.impact == True and args.body_shape == 'Spacecraft'):
    performances_plots(state,m,name,args)
    if args.body_shape == 'Spacecraft' and not cnf.impact:
        closed_form_solution_plot(name,m)
        angle_of_attack_plot(name, args)
    if args.Odyssey_sim == 1:
        ABM_periapsis(name)
            
    plt.close('all')
    # ABM_periapsis(name)

def angle_of_attack_plot(name,args):
    fig , ax = plt.subplots(1,figsize=(21 , 12))
    alt_idx = [idx for idx, val in enumerate(cnf.solution.orientation.alt) if val < args.AE*1e3]

    index_orbit = [0]
    time_0 = [cnf.solution.orientation.time[alt_idx[0]]]
    for i in range(len(alt_idx)):
        if alt_idx[i] - alt_idx[i - 1] > 2:
            index_orbit.append(i)
            time_0.append(cnf.solution.orientation.time[alt_idx[i]])
    index_orbit.append(len(alt_idx))
    if len(index_orbit) == 2: #if 1 orbit
        time = [cnf.solution.orientation.time[i] for i in alt_idx]
        aoa = [np.degrees(cnf.solution.physical_properties.aoa[i]) for i in alt_idx]
        plt.plot(time , aoa , 'k' , linewidth=2)
        plt.xlabel('Time, s' , fontsize=30)
        plt.ylabel('Angle of Attack, degrees' , fontsize=30)
    else:
        time_end = 0
        x_labels = []
        for i in range(0,len(index_orbit) - 1):
            time = [cnf.solution.orientation.time[j] - time_0[i] + time_end for j in
                    alt_idx[index_orbit[i]:index_orbit[i + 1] - 1]]
            aoa = [np.degrees(cnf.solution.physical_properties.aoa[j]) for j in alt_idx[index_orbit[i]:index_orbit[i + 1] - 1]]
            plt.plot(time , aoa , 'k' , linewidth=2)
            plt.plot(time[-1] , aoa[-1] , '.' , color='orangered' , markersize=20)
            time_end = time[-1]
            x_labels.append((time_end+time[0])/2)
        ax.set_xticks(x_labels)
        ax.set_xticklabels(x_labels)

        a = ax.get_xticks().tolist()
        for i in range(0,len(a)):
            a[i] = ' '
        if len(index_orbit)-1<= 5:
            for i in range(0,len(index_orbit)-1):
                a[i] = 'Orbit {}'.format(i+1)
        else:
            s = ((len(index_orbit)-1) + (len(a)-2) // 2) // (len(a)-2)
            interval = (len(a) - 2) // ((len(index_orbit) - 1)//s)
            if interval == 0:
                interval = 1
            r = interval // 2
            for i in range(((len(index_orbit) - 1) + (s) // 2) // (s)):
                step = i*s
                a[i * interval + r ] = '{}'.format(step + 1)
            plt.xlabel('Orbits' , fontsize=30)
        ax.set_xticklabels(a)
        plt.ylabel('Angle of Attack, degrees' , fontsize=30)

    plt.grid()
    plt.tick_params(axis='both' , which='major' , labelsize=20)

    plt.gca().set_ylim([-5 , 100])
    fig.savefig(name + '_angle_of_attack_profile.png')




def closed_form_solution_plot(name,mission):
    fig,(axs1,axs2,axs3) = plt.subplots(1 , 3 , figsize=(22, 10))
    alt = [item - mission.planet.Rp_e for item in cnf.solution.orientation.pos_ii_mag]
    alt_idx = [idx for idx , val in enumerate(alt) if val <= 160*1e3]
    axs1 = plt.subplot(131)
    axs2 = plt.subplot(132)
    axs3 = plt.subplot(133)

    index_orbit = [0]
    time_0 = [cnf.solution.orientation.time[alt_idx[0]]]
    for i in range(len(alt_idx)):
        if alt_idx[i] - alt_idx[i - 1] > 50:
            index_orbit.append(i)
            time_0.append(cnf.solution.orientation.time[alt_idx[i]])
    index_orbit.append(len(alt_idx))


    alt_idx_cf = [idx for idx , val in enumerate(cnf.solution.closed_form.h_cf) if ((val > 0) and (val <= 160*1e3))]
    index_orbit_cf = [0]
    time_0_cf = [cnf.solution.closed_form.t_cf[alt_idx_cf[0]]]
    for i in range(len(alt_idx_cf)):
        if alt_idx_cf[i] - alt_idx_cf[i - 1] > 50:
            index_orbit_cf.append(i)
            time_0_cf.append(cnf.solution.closed_form.t_cf[alt_idx_cf[i]])
    index_orbit_cf.append(len(alt_idx_cf))


    for i in range(0,len(index_orbit)-1):
        time_end = 0
        time = [cnf.solution.orientation.time[j] - time_0[i] + time_end for j in
                alt_idx[index_orbit[i]:index_orbit[i + 1] - 1]]
        alt = [((cnf.solution.orientation.pos_ii_mag[j]-mission.planet.Rp_e)/1e3) for j in
               alt_idx[index_orbit[i]:index_orbit[i + 1] - 1]]
        gamma = [np.degrees(cnf.solution.orientation.gamma_ii[j]) for j in
               alt_idx[index_orbit[i]:index_orbit[i + 1] - 1]]
        v = [(cnf.solution.orientation.vel_ii_mag[j] / 10 ** 3) for j in
               alt_idx[index_orbit[i]:index_orbit[i + 1] - 1]]

        time_cf = [cnf.solution.closed_form.t_cf[j] - time_0_cf[i] + time_end for j in
                   alt_idx_cf[index_orbit_cf[i]:index_orbit_cf[i + 1] - 1]]
        alt_cf = [(cnf.solution.closed_form.h_cf[j]/1e3) for j in
                  alt_idx_cf[index_orbit_cf[i]:index_orbit_cf[i + 1] - 1]]
        gamma_cf = [np.degrees(cnf.solution.closed_form.gamma_cf[j]) for j in
                    alt_idx_cf[index_orbit_cf[i]:index_orbit_cf[i + 1] - 1]]
        v_cf = [(cnf.solution.closed_form.v_cf[j]/1e3) for j in
                alt_idx_cf[index_orbit_cf[i]:index_orbit_cf[i + 1] - 1]]

        closed_form_solution_orbit(time, alt, gamma,v,time_cf,alt_cf,gamma_cf,v_cf, axs1,axs2,axs3)

    axs1.set_xlabel('Time, s', fontsize = 25)
    axs1.set_ylabel('Altitude, km', fontsize = 25)
    axs1.tick_params(labelsize=20)
    axs1.grid()
    plt.subplots_adjust(left=0.08 , right=0.95 , wspace=0.25)

    axs2.set_xlabel('Time, s', fontsize = 25)
    axs2.set_ylabel('Flight-Path Angle, deg', fontsize = 25)
    axs2.tick_params(labelsize=20)
    axs2.grid()

    axs3.set_xlabel('Time, s', fontsize = 25)
    axs3.set_ylabel('Velocity, km/s', fontsize = 25)
    axs3.tick_params(labelsize = 20)
    axs3.grid()

    fig.savefig(name + '_closed_form_solution.png')

def closed_form_solution_orbit(time, alt, gamma,v,time_cf,alt_cf,gamma_cf,v_cf, axs1,axs2,axs3):
    axs1.plot(time,alt, 'k',linewidth=2)
    axs1.plot(time_cf,alt_cf, 'silver',linewidth=2)


    axs2.plot(time,gamma, 'k',linewidth=2)
    axs2.plot(time_cf,gamma_cf , 'silver',linewidth=2)


    axs3.plot(time,v, 'k' , label = 'simulation',linewidth=2)
    axs3.plot(time_cf,v_cf , 'silver',label = 'approx.',linewidth=2)
    axs3.legend(fontsize = 20)


def performances_plots(state, m, name,args):
    fig , axs= plt.subplots(2,2,figsize=(21 , 12))
    if args.body_shape == 'Spacecraft':
        axs[0, 0].plot([item / (60 * 60 * 24) for item in cnf.solution.orientation.time] ,[item / 10 ** 6 for item in cnf.solution.forces.energy] , 'k' , linewidth=2)
        axs[0, 0].set(xlabel='Time, days' , ylabel='Energy, MJ/kg')  # , fontsize=20)
        axs[0, 0].tick_params(axis='both' , which='major' , labelsize=15)
        axs[0, 0].ticklabel_format(useOffset=False)
        axs[0, 0].grid()
        axs[0, 0].xaxis.get_label().set_fontsize(20)
        axs[0, 0].yaxis.get_label().set_fontsize(20)
    else:
        alt = [item for item in cnf.solution.orientation.alt if item < 200*1e3]
        index = [idx for idx, val in enumerate(cnf.solution.orientation.alt) if val < 200*1e3]
        index_orbit = [0]
        time_0 = [cnf.solution.orientation.time[index[0]]]
        for i in range(len(index)):
            if index[i] - index[i - 1] > 50:  # evaluates the number of orbits
                index_orbit.append(i)
                time_0.append(cnf.solution.orientation.time[index[i]])
        index_orbit.append(len(index))
        if len(index_orbit) == 2:
            time = [cnf.solution.orientation.time[i] for i in index]
            axs[0, 0].plot(time, alt, 'k', linewidth=2)
            axs[0, 0].set(xlabel='Time, s', ylabel='Altitude, km')  # , fontsize=20)
        #axs[0, 0].plot([item for item in cnf.solution.orientation.time] , [item/1e3 for item in cnf.solution.orientation.alt],'k',linewidth=2)
        #axs[0, 0].set(xlabel='Time, s', ylabel='Altitude, km')#, fontsize=20)
            axs[0, 0].tick_params(axis='both' , which='major' , labelsize=15)
            axs[0, 0].ticklabel_format(useOffset=False)
            axs[0, 0].grid()
            axs[0, 0].xaxis.get_label().set_fontsize(20)
            axs[0, 0].yaxis.get_label().set_fontsize(20)

    heat_rate = [item for item in cnf.solution.performance.heat_rate if item >0]
    index = [idx for idx, val in enumerate(cnf.solution.performance.heat_rate) if val > 0]
    index_orbit = [0]
    time_0 = [cnf.solution.orientation.time[index[0]]]

    for i in range(len(index)):
        if index[i] - index[i-1]>50:  # evaluates the number of orbits
            index_orbit.append(i)
            time_0.append(cnf.solution.orientation.time[index[i]])
    index_orbit.append(len(index))
    if len(index_orbit) == 2:
        time = [cnf.solution.orientation.time[i] for i in index]
        axs[0 , 1].plot(time , heat_rate , 'k' , linewidth=2)
        axs[0 , 1].set(xlabel='Time, s' , ylabel='Heat Rate, W/cm^2')  # , fontsize=20)
    else:
        time_end = 0
        for i in range(len(index_orbit)-1):
            time = [cnf.solution.orientation.time[j]-time_0[i]+time_end for j in index[index_orbit[i]:index_orbit[i+1]-1]]
            heat_rate = [cnf.solution.performance.heat_rate[j] for j in index[index_orbit[i]:index_orbit[i+1]-1]]
            max_value = max(heat_rate)
            max_index = heat_rate.index(max_value)
#            axs[0,1].plot(time , heat_rate , 'k' , linewidth=2)
            # axs[0,1].plot(time[-1],heat_rate[-1],'.',color='orangered', markersize=12)
            axs[0,1].plot(time[max_index],heat_rate[max_index],'.',color='k', markersize=12)
            time_end = time[-1]
        a = axs[0,1].get_xticks().tolist()
        for i in range(len(a)):
            a[i] = ' '
        if len(index_orbit)-1< (len(a)-2):
            interval = (len(a) - 2) // (len(index_orbit) - 1)
            r = interval // 2
            for i in range(len(index_orbit) - 1):
                a[i * interval + r + 1] = 'Orbit {}'.format(i + 1)
        else:
            s = ((len(index_orbit)-1) + (len(a)-3) // 2) // (len(a)-3)
            interval = (len(a) - 3) // ((len(index_orbit) - 1)//s)
            if interval == 0:
                interval = 1
            r = interval // 2
            for i in range(((len(index_orbit)-1) + (s) // 2) // (s)):
                step = i*s
                a[i * interval + r] = '{}'.format(step + 1)
            axs[0, 1].set(xlabel='Orbits')

        axs[0,1].set_xticklabels(a)
        axs[0,1].set_ylim(bottom=0)
        axs[0,1].set(ylabel='Heat Rate, W/cm^2')  # , fontsize=20)


    axs[0,1].tick_params(axis='both' , which='major' , labelsize=15)
    axs[0,1].grid()
    axs[0,1].xaxis.get_label().set_fontsize(20)
    axs[0,1].yaxis.get_label().set_fontsize(20)

    if len(index_orbit) == 2:
        time = [cnf.solution.orientation.time[i] for i in index]
        heat_load = [cnf.solution.performance.heat_load[i] for i in index]
        axs[1, 0].plot(time , heat_load , 'k' , linewidth=2)
        axs[1 , 0].set(xlabel='Time, s' , ylabel='Heat Load, J/cm^2')
    else:
        time_end = 0
        for i in range(len(index_orbit)-1):
            time = [cnf.solution.orientation.time[j]-time_0[i]+time_end for j in index[index_orbit[i]:index_orbit[i+1]-1]]
            heat_load = [cnf.solution.performance.heat_load[j] for j in index[index_orbit[i]:index_orbit[i+1]-1]]
#            axs[1, 0].plot(time , heat_load , 'k' , linewidth=2)
            axs[1, 0].plot(time[-1] , heat_load[-1] , '.' , color='k' , markersize=12)
            time_end = time[-1]
        a = axs[1,0].get_xticks().tolist()
        for i in range(len(a)):
            a[i] = ' '
        if len(index_orbit)-1< (len(a)-2):
            interval = (len(a) - 2) // (len(index_orbit) - 1)
            r = interval // 2
            for i in range(len(index_orbit) - 1):
                a[i * interval + r + 1] = 'Orbit {}'.format(i + 1)
        else:
            s = ((len(index_orbit)-1) + (len(a)-2) // 2) // (len(a)-2)
            interval = (len(a) - 2) // ((len(index_orbit) - 1)//s)
            if interval == 0:
                interval = 1
            r = interval // 2
            for i in range(((len(index_orbit) - 1) + (s) // 2) // (s)):
                step = i*s
                a[i * interval + r] = '{}'.format(step + 1)
            axs[1 , 0].set(xlabel='Orbits')

        axs[1,0].set_xticklabels(a)
        axs[1 ,0].set(ylabel='Heat Load, J/cm^2')  # , fontsize=20)

    #axs[1,0].plot(time, heat_load,'k',linewidth=2)
    axs[1,0].tick_params(axis='both' , which='major' , labelsize=15)
    axs[1,0].grid()
    axs[1,0].xaxis.get_label().set_fontsize(20)
    axs[1,0].yaxis.get_label().set_fontsize(20)

    if args.body_shape == 'Spacecraft':
        axs[1,1].plot([item / (60 * 60*24) for item in cnf.solution.orientation.time] , [item-args.dry_mass for item in cnf.solution.performance.mass],'k',linewidth=2)
        axs[1,1].set(xlabel='Time, days', ylabel='Fuel Mass, kg')#, fontsize=20)
        axs[1,1].tick_params(axis='both' , which='major' , labelsize=15)
        axs[1,1].ticklabel_format(useOffset=False)
        axs[1,1].grid()
        axs[1,1].xaxis.get_label().set_fontsize(20)
        axs[1,1].yaxis.get_label().set_fontsize(20)
    else:
        # axs[1,1].plot([item for item in cnf.solution.orientation.time] , [item for item in cnf.solution.performance.q],'k',linewidth=2)
        q = [item for item in cnf.solution.performance.q if item > 0.0]
        index = [idx for idx, val in enumerate(cnf.solution.performance.q) if val > 0]
        index_orbit = [0]
        time_0 = [cnf.solution.orientation.time[index[0]]]
        for i in range(len(index)):
            if index[i] - index[i - 1] > 50:  # evaluates the number of orbits
                index_orbit.append(i)
                time_0.append(cnf.solution.orientation.time[index[i]])
        index_orbit.append(len(index))
        if len(index_orbit) == 2:
            time = [cnf.solution.orientation.time[i] for i in index]
            axs[1, 1].plot(time, q, 'k', linewidth=2)
            axs[1, 1].set(xlabel='Time, s', ylabel='Dynamic Pressure, N/m^2')  # , fontsize=20)
            axs[1, 1].tick_params(axis='both', which='major', labelsize=15)
            axs[1, 1].ticklabel_format(useOffset=False)
            axs[1, 1].grid()
            axs[1, 1].xaxis.get_label().set_fontsize(20)
            axs[1, 1].yaxis.get_label().set_fontsize(20)
    fig.savefig(name + '_performances.png')

def traj_2d(state, m, name):
    fig = plt.figure(figsize=(21 , 12))
    circle1 = plt.Circle((0, 0), m.planet.Rp_e/1e3, color='orangered')
    circle2 = plt.Circle((0, 0), m.planet.Rp_e/1e3+160, color='y', alpha=0.25)
    x = [item/1e3 for item in cnf.solution.orientation.pos_ii[0]]
    y = [item/1e3 for item in cnf.solution.orientation.pos_ii[1]]
    z =  [item/1e3 for item in cnf.solution.orientation.pos_ii[2]]
    i = cnf.solution.orientation.oe[2][0]
    OMEGA = cnf.solution.orientation.oe[3][0]
    omega = cnf.solution.orientation.oe[4][0]
    T_ijk = [[np.cos(OMEGA) * np.cos(omega) - np.sin(OMEGA) * np.sin(omega) * np.cos(i),
             np.sin(OMEGA) * np.cos(omega) + np.cos(OMEGA) * np.sin(omega) * np.cos(i),
             np.sin(omega) * np.sin(i)],
             [-np.cos(OMEGA) * np.sin(omega) - np.sin(OMEGA) * np.cos(omega) * np.cos(i),
              -np.sin(OMEGA) * np.sin(omega) + np.cos(OMEGA) * np.cos(omega) * np.cos(i),
              np.cos(omega) * np.sin(i)],
             [np.sin(OMEGA) * np.sin(i), -np.cos(OMEGA) * np.sin(i), np.cos(i)]]
    vector = np.inner(T_ijk,np.c_[x,y,z])

    plt.plot(vector[0],vector[1],'k',linewidth=2)
    plt.plot(vector[0][0],vector[1][0],'.',c='red')
    ax = fig.gca()
    ax.add_artist(circle1)
    ax.add_artist(circle2)
    ax.set_xlabel('X, km', fontsize=20, labelpad=15)
    ax.set_xlim(-state['Apoapsis']/1e3, state['Periapsis']/1e3)
    ax.set_ylabel('Y, km', fontsize=20, labelpad=15)
    ax.set_ylim(-3*m.planet.Rp_e/1e3, 3*m.planet.Rp_e/1e3)
    ax.tick_params(axis='both' , which='major' , labelsize=15)
    fig.savefig(name + '_traj2d.png')

def traj_3d(state,m, name):
    fig = plt.figure(figsize=(21 , 12))

    gs = GridSpec(3, 3 , left=0.05, right=0.98, wspace=0.07)
    ax1 = fig.add_subplot(gs[:, :-1],projection='3d')
    ax2 = fig.add_subplot(gs[-3 , -1],projection='3d')
    ax3 = fig.add_subplot(gs[-2 , -1],projection='3d')
    ax4 = fig.add_subplot(gs[-1 , -1],projection='3d')

    x = [item/1e3 for item in cnf.solution.orientation.pos_ii[0]]
    y = [item/1e3 for item in cnf.solution.orientation.pos_ii[1]]
    z = [item/1e3 for item in cnf.solution.orientation.pos_ii[2]]

    stride = 2

    u = np.linspace(0 , 2 * np.pi , 100)
    v = np.linspace(0 , np.pi , 100)
    r = m.planet.Rp_e / 10 ** 3
    x_s = r * np.outer(np.cos(u) , np.sin(v))
    y_s = r * np.outer(np.sin(u) , np.sin(v))
    z_s = r * np.outer(np.ones(np.size(u)) , np.cos(v))
    ax1.plot_surface(x_s , y_s , z_s , linewidth=0.0 , cstride=stride , rstride=stride , color='orangered' ,
                     edgecolors='k' , lw=0.6 , zorder=2)
    ax2.plot_surface(x_s , y_s , z_s , linewidth=0.0 , cstride=stride , rstride=stride , color='orangered' ,
                     edgecolors='k' , lw=0.6 , zorder=2)
    ax3.plot_surface(x_s , y_s , z_s , linewidth=0.0 , cstride=stride , rstride=stride , color='orangered' ,
                     edgecolors='k' , lw=0.6 , zorder=2)
    ax4.plot_surface(x_s , y_s , z_s , linewidth=0.0 , cstride=stride , rstride=stride , color='orangered' ,
                     edgecolors='k' , lw=0.6 , zorder=2)

    r += 160
    x_s = r * np.outer(np.cos(u) , np.sin(v))
    y_s = r * np.outer(np.sin(u) , np.sin(v))
    z_s = r * np.outer(np.ones(np.size(u)) , np.cos(v))
    ax1.plot_surface(x_s , y_s , z_s , rstride=4 , cstride=4 , color='y' , linewidth=0 , alpha=0.25 , zorder=2)
    ax2.plot_surface(x_s , y_s , z_s , rstride=4 , cstride=4 , color='y' , linewidth=0 , alpha=0.25 , zorder=2)
    ax3.plot_surface(x_s , y_s , z_s , rstride=4 , cstride=4 , color='y' , linewidth=0 , alpha=0.25 , zorder=2)
    ax4.plot_surface(x_s , y_s , z_s , rstride=4 , cstride=4 , color='y' , linewidth=0 , alpha=0.25 , zorder=2)

    index = [0]
    for i in range(len(cnf.solution.orientation.numberofpassage)-1):
        if (cnf.solution.orientation.numberofpassage[i+1]-cnf.solution.orientation.numberofpassage[i]>0):
            index.append(i)
    index.append(len(cnf.solution.orientation.numberofpassage))

    for i in range(len(index)-1):
        xs = x[index[i]:index[i+1]]
        ys = y[index[i]:index[i+1]]
        zs = z[index[i]:index[i+1]]
        traj_3d_orbit(xs,ys,zs,state,ax1,ax2,ax3,ax4)

    fig.subplots_adjust(left=0 , right=1 , bottom=0 , top=1)
    fig.savefig(name + '_traj3d.png')
    return

def traj_3d_orbit(x,y,z,state,ax1,ax2,ax3,ax4):
    ax1.plot(x , y , z , 'k' , linewidth=1,zorder=1)
    ax2.plot(x , y , z , 'k' , linewidth=1,zorder=1)
    ax3.plot(x , y , z , 'k' , linewidth=1,zorder=1)
    ax4.plot(x , y , z , 'k' , linewidth=1,zorder=1)



    a = (state['Apoapsis'] / 10 ** 3 + state['Periapsis'] / 10 ** 3 + 3396) / 2
    e = (state['Apoapsis'] / 10 ** 3 - state['Periapsis'] / 10 ** 3 - 3396) / (
                state['Apoapsis'] / 10 ** 3 + state['Periapsis'] / 10 ** 3 + 3396)
    b = a * (1 - e ** 2) ** 0.5
    c = (+state['Apoapsis'] / 10 ** 3 + state['Periapsis'] / 10 ** 3)/2
    ax1.set_xlabel('X, km' , fontsize=15 , labelpad=15)
    ax1.set_xlim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    ax1.set_ylabel('Y, km' , fontsize=15 , labelpad=15)
    ax1.set_ylim(-c,c)
    # ax1.set_ylim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    ax1.set_zlabel('Z, km' , fontsize=15 , labelpad=15)
    ax1.set_zlim(-c,c)
    # ax1.set_zlim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    ax1.view_init(220 , 225)
    x1 , y1 , z1 = plot_visible(220 , 225 , x , y , z)
    plots_continuous(x1 , y1 , z1 , ax1)
    ax1.tick_params(axis='both' , which='major' , labelsize=12)
    ax1.dist = 10

    # ax2.set_xlabel('X, km', fontsize=15)
    ax2.set_xlim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    plt.setp(ax2.get_xticklabels() , visible=False)
    ax2.set_ylabel('Y, km' , fontsize=15 , labelpad=15)
    ax2.set_ylim(-c,c)
    # ax2.set_ylim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    ax2.set_zlabel('Z, km' , fontsize=15 , labelpad=15)
    # ax2.set_zlim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    ax2.set_zlim(-c,c)
    ax2.view_init(0 , 0)
    x2 , y2 , z2 = plot_visible(0 , 0 , x , y , z)
    plots_continuous(x2 , y2 , z2 , ax2)
    ax2.dist = 8

    ax3.set_xlabel('X, km' , fontsize=15 , labelpad=15)
    ax3.set_xlim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    # ax3.set_ylabel('Y, km', fontsize=15)
    # ax3.set_ylim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    ax3.set_ylim(-c,c)
    plt.setp(ax3.get_yticklabels() , visible=False)
    ax3.set_zlabel('Z, km' , fontsize=15 , labelpad=15)
    # ax3.set_zlim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    ax3.set_zlim(-c,c)
    x3 , y3 , z3 = plot_visible(0 , -90 , x , y , z)
    plots_continuous(x3 , y3 , z3 , ax3)
    ax3.view_init(0 , -90)
    ax3.dist = 8

    ax4.set_xlabel('X, km' , fontsize=15 , labelpad=15)
    ax4.set_xlim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    ax4.set_ylabel('Y, km' , fontsize=15 , labelpad=15)
    ax4.set_ylim(-c,c)
    # ax4.set_ylim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    # ax4.set_zlabel('Z, km', fontsize=15)
    # ax4.set_zlim(-state['Apoapsis'] / 10 ** 3 , state['Periapsis'] / 10 ** 3)
    ax4.set_zlim(-c,c)
    plt.setp(ax4.get_zticklabels() , visible=False)
    ax4.view_init(-90 , 0)
    x4 , y4 , z4 = plot_visible(-90 , 0 , x , y , z)
    plots_continuous(x4 , y4 , z4 , ax4)
    ax4.dist = 8

def plot_visible(elev, azimuth,x,y,z):
    #transform viewing angle to normal vector in data coordinates
    a = azimuth*np.pi/180. -np.pi
    e = elev*np.pi/180. - np.pi/2.
    X = [ np.sin(e) * np.cos(a),np.sin(e) * np.sin(a),np.cos(e)]
    # concatenate coordinates
    Z = np.c_[x,y,z]

    # calculate dot product
    # the points where this is positive are to be shown

    cond =(np.dot(Z,X) >= 0)# np.logical_or(np.dot(Z1,X) >= 0, -np.linalg.norm(Z1,axis=1)*np.cos(e)*np.cos(a)>=radius_sphere)
    # filter points by the above condition
    x_ca = np.array(x)[cond]
    y_ca = np.array(y)[cond]
    z_ca = np.array(z)[cond]

    return x_ca, y_ca, z_ca

def plots_continuous(x,y,z,ax):
    index_disc = [0]
    for i in range(1, len(x)):
        if ((x[i] - x[i - 1]) ** 2 +(y[i] - y[i - 1])**2 + (z[i] - z[i - 1])**2)**0.5> 500:
            index_disc.append(i-1)
            index_disc.append(i)
    index_disc.append(len(x))
    for i in range(0, len(index_disc)-1, 2):
        ax.plot(x[index_disc[i]:index_disc[i+1]] , y[index_disc[i]:index_disc[i+1]], z[index_disc[i]:index_disc[i+1]] , 'k' ,  linewidth=1, zorder=5)


def ABM_periapsis(name):
    # Create some mock data
    orbit_number = cnf.orbit_number_list
    periapis_altitude = cnf.periapsis_list
    delta_v = [0]
    raise_man_orbit = [7,55,87,110,127,160,178,193,210,222]
    delta_v_raise = []
    lower_man_orbit = [13,25,30,36,48,68,72,80,]
    delta_v_lower = []
    for j in range(0,len(cnf.delta_v_list)-1):
        delta_v.append(cnf.delta_v_list[j]-cnf.delta_v_list[j-1])

    for j in range(len(delta_v)):
        if j in raise_man_orbit:
            delta_v_raise.append(delta_v[j])
        elif j in lower_man_orbit:
            delta_v_lower.append(delta_v[j])
    fig, ax1 = plt.subplots()


    ax1.set_xlabel('Orbit Number', fontsize=14)
    ax1.set_ylabel('Periapsis Altitude, km', color='k', fontsize=14)
    ax1.plot(orbit_number , periapis_altitude , '.', color='k', markersize = 12, label = 'Periapsis Altitude')
    ax1.tick_params(axis='both', labelcolor='k', labelsize=12)
    ax1.set_ylim([90,160])
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis



    ax2.set_ylabel('ABM magnitude, m/s' , color='k',fontsize = 14)  # we already handled the x-label with ax1
    ax2.plot([item +19 for item in  raise_man_orbit], delta_v_raise ,"^", color='k', markersize = 12, label='ABM to raise periapsis')
    ax2.plot([item +19 for item in lower_man_orbit], delta_v_lower ,"o",color='k', markersize = 12, label='ABM to lower periapsis')
    ax2.tick_params(axis='y', labelcolor='k', labelsize=12)
    ax2.set_ylim([0.0,1.4])
    ax1.grid()

    # ask matplotlib for the plotted objects and their labels
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=0, fontsize = 'large')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    fig.savefig(name + '_odyssey.png')
