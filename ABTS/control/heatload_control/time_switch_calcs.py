from scipy import optimize
def switch_calculation_with_integration(ip, m,position,args,t,heat_rate_control,reevaluation_mode, current_position = 0):
    from physical_models.Aerodynamic_models import aerodynamicscoefficient_fM as am
    import scipy
    import numpy as np
    import time
    from control.eoms import asim
    # start_time = time.time()

    def func(k_cf):
        global time_switch, k
        k = k_cf/100
        [y , time_switch] = asim(ip, m , t , position , args , k  , heat_rate_control, time_switch_eval = True)
        # print('k_cf',k)
        Q = y[-1][-1]
        # print('Delta Q',(Q - m.aerodynamics.heat_load_limit))
        return (Q - m.aerodynamics.heat_load_limit)


    # Evaluates max Q
    delta_Q_max = func(5)
    delta_Q_min = func(0)

    # print('delta_Q_max',delta_Q_max,'delta_Q_min',delta_Q_min)
    if delta_Q_max*delta_Q_min <0.0:
        # try:
        optimize.bisect(func,3.3, 4.5, xtol =1e-7)
        # except: # if on the boundary and do not find intersections
        #     print('here')
        #     optimize.bisect(func, k*100-0.01, k*100+0.01, xtol =1e-2)
    elif delta_Q_max < 0.0:
        if (args.heat_load_sol == 0 or args.heat_load_sol == 3):
            return [0.0,0.0]
        elif (args.heat_load_sol == 1 or args.heat_load_sol == 2):
            return [0.0, 1000.0]
    elif delta_Q_min > 0.0:
        if (args.heat_load_sol == 0 or args.heat_load_sol == 3):
            return [0.0, 1000.0]
        elif (args.heat_load_sol == 1 or args.heat_load_sol == 2):
            return [0.0, 0.0]
    # print("--- %s seconds ---" % (time.time() - start_time))
    return time_switch

def switch_calculation(ip, m,position,args,t,heat_rate_control,reevaluation_mode,current_position = 0):
    from physical_models.Aerodynamic_models import aerodynamicscoefficient_fM as am
    import scipy
    from control.Control import heat_rate_calc
    from control.heatload_control.utils_timeswitch import lambdas
    from control.heatload_control.utils_timeswitch import aoa
    from control.heatload_control.utils_timeswitch import func
    from utils.closed_form_solution import closed_form
    import math
    from physical_models.Density_models import density_exp as dm
    import config
    import numpy as np


    global t_cf, h_cf, gamma_cf, v_cf



    # Evaluate initial conditions
    T = m.planet.T  # fixed temperature
    t_cf, h_cf, gamma_cf, v_cf = closed_form(args, m, position, T, aoa=m.aerodynamics.aoa,
                                                 online=True)  # define closed-form solution
    RT = T * m.planet.R
    S = (v_cf / (2 * RT) ** 0.5)
    CL_90, CD_90 = am(np.pi / 2, m.body, T, S[0], m.aerodynamics,
                          montecarlo=0)  # Cl and Cd for all the body perpendicular
    CL_0, CD_0 = am(0, m.body, T, S[0], m.aerodynamics, montecarlo=0)  # Cl and Cd for all the body perpendicular
    CD_slope = (CD_90 - CD_0) / (math.pi * 0.5)
    coeff = [CD_slope,CL_0, CD_0]


    # Evaluates max Q
    aoa_cf = aoa(m,0.1, t_cf , h_cf , gamma_cf , v_cf,coeff)[0]
    approx_sol = [t_cf , h_cf , gamma_cf , v_cf]
    delta_Q_max = func(0.1, m, args, coeff, position, heat_rate_control, approx_sol,   aoa_cf)
    delta_Q_min = func(0.0, m, args, coeff, position, heat_rate_control, approx_sol, np.array([0.0]*len(aoa_cf)))
    if delta_Q_max*delta_Q_min <0.0:
        k_cf = optimize.brentq(func ,0.0 , 0.1,xtol=1e-7, args =(m,args, coeff,position,heat_rate_control,approx_sol,aoa_cf)) #,xtol=1e-6,
    elif delta_Q_max < 0.0:
        return [0.0,0.0]
    elif delta_Q_min > 0:
        return [0.0, t_cf[-1]/2]#t_cf[-1]-100]

    [t_cf,v_cf,gamma_cf,h_cf] = func(k_cf,m,args,coeff, position, heat_rate_control,approx_sol, aoa_cf, approx_calc=True)

    lambda_switch, lambdav = lambdas(m,aoa_cf, k_cf, t_cf , h_cf , gamma_cf , v_cf,coeff)[0:2]

    index_array = (lambdav < lambda_switch)
    temp = t_cf*index_array
    temp = np.trim_zeros(temp)
    t_switch = [temp[0],temp[-1]]
    if abs(t_switch[-1] - t_cf[-1])<5:
        t_switch[0] -= t_cf[-1]*0.04
    elif abs(t_switch[-1] - t_cf[-1])<60:
        t_switch[0] -= 5
    t_switch[1] -= t_switch[1]*0.1
    return t_switch