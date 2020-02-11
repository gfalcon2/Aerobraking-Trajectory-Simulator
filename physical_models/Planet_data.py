from config import model as m
def planet_data(ip):
    try:
        if 'Planet' in ip.keys():
            ip = ip['Planet']
    except:
        pass

    if (ip == 0): # Earth
        Rp_e =  6.3781e6  # m, equatorial radius
        Rp_p = 6.3568e6  # polar radius, m
        Rp_m = 6.3710e6  # volumetric mean radius, m
        mass = 5.9736e24  # mass, kg
        g_ref = 9.798  # m/s^2
        rho_ref = 1.225  # kg/m^3
        mu = 3.9860e14  # gravitational parameter, m^3/s^2
        h_ref = 0 * 10 ** 3
        H = 8.5 * 10 ** 3  # m
        R = 287.1  # J/KgK
        gamma = 1.4005  # ratio of specific heats, nd
        T = 300 # K constant
        p = 101400  # kg/ms, Surface pressure
        J2 = 1.0826e-3
        k = 1.83e-4  # Chapman heating coefficient, kg ^ 0.5 / m
        # k = 1.7623e-4 # Sutton - Graves heating coefficient, kg^0.5 / m
        omega = [0, 0,  7.2921066e-5]  # angular velocity vector, rad / s
        mu_fluid = 1.5*10e-5#m2 s−1 Kinematic viscosity
        Lz = -9.8/10**3 # K/m


    if (ip == 1):  # Mars
        Rp_e = 3.3962e6  # m
        Rp_p = 3.3762e6 # polar radius, m
        Rp_m = 3.3895e6  # volumetric mean radius, m
        mass = 6.4185e23  # mass, kg
        g_ref = 3.71  # m/s^2
        rho_ref = 3.8*10**-8#7.3*10**-8#0.02  # kg/m^3
        mu = 4.2828e13 # gravitational parameter, m^3/s^2
        h_ref = 103*10**3#109 * 10 ** 3
        H =9 * 10 ** 3 #10.6 * 10 ** 3  # m
        R = 188.92  # J/KgK  # wrong check
        gamma = 1.33  # wrong check
        T = 190  # K constant # wrong anchor # the script is not considering this number but T =150
        p = 636  # N/m^2, Surface pressure # wrong anchor
        J2 = 1.96045e-3
        k = 1.898e-4  # Sutton - Graves heating coefficient, kg ^ 0.5 / m
        omega = [0, 0,  7.0882360e-5]  # angular velocity vector, rad / s
        mu_fluid = 13.06*10e-6 # Pa*s  Kinematic viscosity
        Lz = -4.5/10**3 # K/m

    if (ip == 2):  # Venus
        Rp_e = 6.0518e6 # equatorial radius, m
        Rp_p = 6.0518e6 # polar radius, m
        Rp_m = 6.0518e6 # volumetric mean radius, m
        mass = 4.8685e24 # mass, kg
        g_ref = 8.87  # m/s^2
        rho_ref = 65  # kg/m^3
        mu = 3.249e14 # gravitational parameter, m^3/s^2
        h_ref = 0 * 10 ** 3
        H = 15.9 * 10 ** 3  # m
        R = 188.92  # J/KgK  # ideal gas
        gamma = 1.2857 # ratio of specific heats, nd
        T = 100  # K, reference temperature constant
        p = 9200000 # kg/ms, Surface pressure
        J2 = 4.458e-6
        k = 1.896e-4 # Sutton - Graves heating coefficient, kg ^ 0.5 / m
        omega = [0, 0, -2.9924205e-7] # angular velocity vector, rad / s
        Lz = -10.7/10**3 # K/m



    return m.planet(Rp_e, Rp_p, Rp_m, mass, p, k, omega, g_ref, rho_ref, h_ref, H, R, gamma, T, J2, mu, mu_fluid, Lz)