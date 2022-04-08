import math
import numpy as np
from physical_models.Density_models import density_exp as dm
from control.Control import heat_rate_calc
from utils.closed_form_solution import closed_form
def lambdas(m,aoa, k, t, h, gamma, v, coeff):
    CD_slope, CL_0, CD_0 = coeff
    S_ref = m.body.Area_tot
    mass = m.body.Mass
    Rp = m.planet.Rp_e
    mu = m.planet.mu
    g0 = m.planet.g_ref
    H = m.planet.H
    # evaluate lambda_switch
    lambda_switch = np.divide(k * 2.0 * mass * v,
                              S_ref * CD_slope * math.pi)
    lambdav = [0] * len(t)
    lambdag = [0] * len(t)
    lambdah = [0] * len(t)
    lambdav[-1] = v[-1]
    lambdag[-1] = 0.0
    lambdah[-1] = mu / (Rp + h[-1]) ** 2
    rho = dm(h=h, p=m.planet)[0]
    # g = np.divide(mu,np.square(Rp+h))
    g = np.divide(g0 * Rp ** 2, np.square(Rp + h))  # np.divide(-mu,np.square(Rp+h))
    for ii in range(len(t) - 1, 0, -1):  # reversed(t):
        # print(ii)
        lambdav_dot = -3 * k * rho[ii - 1] * v[ii - 1] ** 2 * aoa[ii - 1] / math.pi + lambdav[ii] * (
                    rho[ii - 1] * S_ref * (CD_0 + aoa[ii - 1] * CD_slope) * v[ii - 1]) / mass - \
                      lambdag[ii] * ((rho[ii - 1] * S_ref * CL_0) / (2 * mass) + g[ii - 1] / v[ii - 1] ** 2 + 1 / (
                    Rp + h[ii - 1])) - lambdah[ii] * gamma[ii - 1]
        lambdag_dot = lambdav[ii] * g[ii - 1] - lambdah[ii] * v[ii - 1]
        lambdah_dot = k * rho[ii - 1] * v[ii - 1] ** 3 * aoa[ii - 1] / (math.pi * H) - lambdav[ii] * (
                    (rho[ii - 1] * S_ref * (CD_0 + aoa[ii - 1] * CD_slope) * v[ii - 1] ** 2) / (2 * mass * H) + 2 * g[
                ii - 1] * gamma[ii - 1] / (Rp + h[ii - 1])) \
                      + lambdag[ii] * (rho[ii - 1] * S_ref * CL_0 * v[ii - 1] / (2 * mass * H) - 2 * g[ii - 1] / (
                    (Rp + h[ii - 1]) * v[ii - 1]) + v[ii - 1] / (Rp + h[ii - 1]) ** 2)
        lambdav[ii - 1] = lambdav[ii] - lambdav_dot * (t[ii] - t[ii - 1])
        # print(lambdav[ii-1])
        lambdag[ii - 1] = lambdag[ii] - lambdag_dot * (t[ii] - t[ii - 1])
        lambdah[ii - 1] = lambdah[ii] - lambdah_dot * (t[ii] - t[ii - 1])

    in_cond_lambda = [lambdav[1], lambdag[1], lambdah[1]]
    # import matplotlib.pyplot as plt
    #
    # plt.plot(t,lambdav)
    # plt.plot(t,lambda_switch)
    # plt.show()
    return lambda_switch, lambdav, in_cond_lambda

def aoa(m, k_cf, t_cf , h_cf , gamma_cf , v_cf, coeff, aoa_cf=[]):
    if not aoa_cf:
        aoa_cf = [m.aerodynamics.aoa] * len(t_cf)

    lambda_switch , lambda_v,in_cond_lambda = lambdas(m,aoa_cf, k_cf,t_cf , h_cf , gamma_cf , v_cf,coeff)

    aoa_cf = [m.aerodynamics.aoa]*len(lambda_v)
    index_array = (lambda_v >= lambda_switch)
    aoa_cf = aoa_cf*index_array
    return np.asarray(aoa_cf), in_cond_lambda

def func(k_cf,m,args,coeff, position, heat_rate_control,approx_sol, aoa_cf, initial_guess=False, approx_calc=False):
    t_cf , h_cf , gamma_cf , v_cf = approx_sol
    # global t_cf , h_cf , gamma_cf , v_cf
    # initialize heat load
    temp_Q = 1000
    Q_prec = 0

    count = 0
    while temp_Q > 1e-3: # stop if two following trajectory evaluation provides the same Q
    # define angle of attack lagrangian multipliers
        T = m.planet.T  # fixed temperature
        [aoa_cf,in_cond_lambda] = aoa(m,k_cf, t_cf , h_cf , gamma_cf , v_cf, coeff, np.ndarray.tolist(aoa_cf))  # update angle of attack profile with new k
        t_cf , h_cf , gamma_cf , v_cf = closed_form(args , m , position , T , aoa=m.aerodynamics.aoa , online=True ,
                                                aoa_profile=np.ndarray.tolist(aoa_cf))  # re-evaluate the closed form solution using previous angle of attack profile
        a = (m.planet.gamma * m.planet.R * T) ** 0.5
        M = v_cf / a
        S = (m.planet.gamma*0.5) ** 0.5 * M
        rho = dm(h=h_cf, p=m.planet)[0]  # density calculated through exponential density

        heat_rate = heat_rate_calc(args.multiplicative_factor_heatload * m.aerodynamics.thermal_accomodation_factor ,
                               rho , T , T , m.planet.R , m.planet.gamma , S ,
                               aoa_cf)  # (list(map(heat_rate_calc, m.aerodynamics.thermal_accomodation_factor,rho,T,T,m.planet.R,m.planet.gamma, S, aoa_cf)))

        if heat_rate_control == True:# Account for max heat rate possible
            index_hr = heat_rate > args.max_heat_rate
            heat_rate[index_hr] =  args.max_heat_rate

        #heat_rate[0:len(config.heat_rate_list)] = config.heat_rate_list
        Q = sum(heat_rate) * t_cf[-1] / len(t_cf)

        # update error
        temp_Q = Q-Q_prec
        # print('temp_Q', temp_Q, 'Q_prec', Q_prec,'Q',Q, 'count', count)
        Q_prec = Q
        count +=1

    if initial_guess:
        return in_cond_lambda
    elif approx_calc:
        return t_cf,v_cf,gamma_cf,h_cf

    delta_Q = (Q - m.aerodynamics.heat_load_limit)
    # print(delta_Q)
    return delta_Q