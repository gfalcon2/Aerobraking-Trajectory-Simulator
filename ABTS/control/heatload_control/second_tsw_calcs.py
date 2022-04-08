import config
import math
from utils.closed_form_solution import closed_form
import numpy as np
from scipy import optimize

def second_time_switch_recalc_with_integration(ip, m,position,args,t, heat_rate_control,reevaluation_mode,current_position = 0):
    # print("Begin second time switch recalculation")
    time_switch = config.time_switch_2
    # Q_past = sum(heat_rate_past)/args.trajectory_rate
    from control.eoms import asim
    def func(t_s):
        # print('t_s',t_s)
        y = asim(ip, m , t , current_position , args, 0, heat_rate_control, t_s,reevaluation_mode)

        Q = y[-1][-1]
        # print('t_s',t_s,'Q in func',Q)
        return (Q - m.aerodynamics.heat_load_limit)

    # Evaluates max Q
    delta_Q_switch = func(time_switch)
    delta_Q_current_time = func(t)


    delta_Q_with_current_time = abs(delta_Q_current_time-delta_Q_switch)# if the difference between the heat load at the current time and the switch is small but the time is too large, recheck

    if abs(delta_Q_switch) <= 0.01: #0.1
            if (args.heat_load_sol == 0 or args.heat_load_sol == 3) and not (delta_Q_with_current_time<0.5 and abs(t-time_switch)>20) and (delta_Q_switch<delta_Q_current_time):
                return  config.time_switch_1, config.time_switch_2
            elif (args.heat_load_sol == 1 or args.heat_load_sol == 2):
                return  config.time_switch_1, config.time_switch_2
    elif delta_Q_current_time < 0:
        if args.heat_load_sol == 0 or args.heat_load_sol == 3:
            return  config.time_switch_1, t
        # if args.heat_load_sol == 2:
        #     return [0.0, 1000.0]

    if reevaluation_mode == 1:
        x_tol = 0.01 #0.1
    else:
        x_tol = 0.01

    if args.heat_load_sol == 0 or args.heat_load_sol == 1:
        b = t +  200
    else:
        b = t+ 1500
    try:
        # print('begin brentq between',t,'and',b)
        time_switch = optimize.brentq(func, t, b, xtol=x_tol)
        # print('new time switch',time_switch)
    except:
        pass
    # print("End second time switch recalculation")

    return config.time_switch_1,time_switch

def second_time_switch_recalc(ip, m,position,args,t, heat_rate_control, current_position = 0, reevaluation_mode = 0):
    # print("Begin second time switch recalculation")
    from physical_models.Density_models import density_exp as dm
    from control.Control import heat_rate_calc
    # print("Begin second time switch recalculation")

    # evaluates past heat load
    aoa_past = config.aoa_list
    time_switch_1 = config.time_switch_1
    time_switch_2 = config.time_switch_2
    # print('t',t,'time switch 2',time_switch_2)
    Q_past = config.heat_load_past

    def func(time_switch):
        # predict the rest part of the passage heat load
        T = m.planet.T #fixed temperature
        t_cf =closed_form(args,m,position,T,aoa=m.aerodynamics.aoa,online=True)[0] # closed-form solution only for length of t

        mask = (t_cf < time_switch_1) | (t_cf > time_switch)
        aoa_list = m.aerodynamics.aoa*mask
        aoa_list = aoa_list.tolist()

        t_cf, h_cf, gamma_cf, v_cf =closed_form(args,m,position,T,aoa_profile=aoa_list,online=True) # closed form solution with the real aoa (this is real only if control mode 2, if control mode 3, angle of attack taken as pi/2 instead of the correct one)
        RT = T*m.planet.R
        S = (v_cf/(2*RT)**0.5)
        index_tilltsw = (t_cf > t) & (t_cf <= time_switch) #[index for index in range(len(t_cf)) if (t_cf[index]<=time_switch_2 and t_cf[index]>t)]
        index_remaining = t_cf > time_switch # [index for index in range(len(t_cf)) if (t_cf[index]>time_switch_2)]
        t_tilltsw = t_cf * index_tilltsw
        t_tilltsw = np.trim_zeros(t_tilltsw)
        t_remaining =  t_cf * index_remaining
        t_remaining = np.trim_zeros(t_remaining)

        if len(t_tilltsw)+len(t_remaining)==0:
            return Q_past
        rho = dm(h=h_cf, p=m.planet)[0]
        rho_tilltsw = rho[index_tilltsw]
        rho_remaining = rho[index_remaining]

        aoa_cf = [0]*len(t_tilltsw)
        S_tilltsw = S[index_tilltsw]
        heat_rate_tilltsw = heat_rate_calc(args.multiplicative_factor_heatload * m.aerodynamics.thermal_accomodation_factor , rho_tilltsw , T , T , m.planet.R , m.planet.gamma , S_tilltsw , aoa_cf)
        if len(t_tilltsw)>1:
            Q_rate_tilltsw=sum(heat_rate_tilltsw)*(t_tilltsw[-1]-t_tilltsw[0])/len(t_tilltsw)
        else:
            Q_rate_tilltsw = 0

        S_remaining = S[index_remaining]
        aoa_cf = [m.aerodynamics.aoa] * len(t_remaining)
        heat_rate_remaining = heat_rate_calc(args.multiplicative_factor_heatload * m.aerodynamics.thermal_accomodation_factor , rho_remaining , T , T , m.planet.R , m.planet.gamma , S_remaining , aoa_cf)

        if args.control_mode == 3: # Account for max heat rate possible
            index_hr = heat_rate_remaining > args.max_heat_rate
            index_hr_not = heat_rate_remaining <= args.max_heat_rate
            temp = args.max_heat_rate*index_hr
            heat_rate_remaining= heat_rate_remaining*index_hr_not+temp
            heat_rate_remaining = heat_rate_remaining.tolist()

        if len(t_remaining) > 1:
            Q_rate_remaining = sum(heat_rate_remaining) * (t_remaining[-1] - t_remaining[0]) / len(t_remaining)
        else:
            Q_rate_remaining = 0

        Q = Q_past+Q_rate_remaining+Q_rate_tilltsw
        return (Q-m.aerodynamics.heat_load_limit)

    # Evaluates max Q
    delta_Q_switch = func(time_switch_2)
    delta_Q_current_time = func(t)
    if reevaluation_mode == 1:
        x_tol = 0.1
    else:
        x_tol = 0.01
    # print('delta_Q_current_time',delta_Q_current_time, 'delta_Q_switch',delta_Q_switch)

    delta_Q_with_current_time = abs(delta_Q_current_time-delta_Q_switch)# if the difference between the heat load at the current time and the switch is small but the time is too large, recheck
    if abs(delta_Q_switch) <= x_tol and not (delta_Q_with_current_time<0.5 and abs(t-time_switch_2)>20):
            return  config.time_switch_1,config.time_switch_2
    elif delta_Q_current_time < 0:
        return  config.time_switch_1,t
    elif delta_Q_switch > 0:
        mult = 1
    elif delta_Q_switch < 0:
        mult = -1

    if reevaluation_mode == 1:
        x_tol = 0.1
    else:
        x_tol = 0.05

    b = t +  200
    try:
        # print('begin brentq between',t,'and',b)
        time_switch_2 = optimize.brentq(func, t, b, xtol=x_tol)
        # print('new time switch',time_switch_2)
    except:
        pass
    time_switch_2 -= time_switch_2*0.1
    return config.time_switch_1,time_switch_2
