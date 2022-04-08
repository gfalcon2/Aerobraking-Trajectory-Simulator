def security_mode(ip,m,position,args,t,heat_rate_control=False):
    from physical_models.Aerodynamic_models import aerodynamicscoefficient_fM as am
    import scipy
    import config
    from control.Control import heat_rate_calc
    from utils.closed_form_solution import closed_form
    import math
    from physical_models.Density_models import density_exp as dm

    T = m.planet.T #fixed temperature
    t_cf, h_cf, gamma_cf, v_cf =closed_form(args,m,position,T,aoa=m.aerodynamics.aoa,online=True) # define closed-form solution
    RT = T*m.planet.R
    S = (v_cf/(2*RT)**0.5)
    rho = dm(h=h_cf, p=m.planet)[0]  # density calculated through exponential density
    # Security Mode
    aoa_cf_min = [0]*len(t_cf)
    heat_rate_min = heat_rate_calc(args.multiplicative_factor_heatload * m.aerodynamics.thermal_accomodation_factor ,
                                   rho ,
                                    T , T , m.planet.R , m.planet.gamma , S ,
                                       aoa_cf_min)

    index_future = t_cf > t
    heat_rate_min = heat_rate_min*index_future
    traj_rate = (t_cf[1]-t_cf[0])
    # print(t,config.heat_load_past,(sum(heat_rate_min.tolist()) * traj_rate)+config.heat_load_past)
    if (sum(heat_rate_min.tolist()) * traj_rate)+config.heat_load_past>m.aerodynamics.heat_load_limit:
        # print('here')
        config.security_mode = True
        # exit()
        return [0, t_cf[-1]+10000.0] # switch angle of attack to 0 for the rest # security made on the check

    return [config.time_switch_1,config.time_switch_2]
