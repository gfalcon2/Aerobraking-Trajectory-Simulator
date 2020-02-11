import math


def ic_calculation_rptoae(planet, gamma, v):
    ra_list = []
    hp_list = []
    for v0 in v:
        for gamma0 in gamma:
            r = planet.Rp_e + 160*10**3  # Drag passage always start and end at 160 km of altitude
            a = planet.mu /((2*planet.mu/r)-v0**2)
            h = r*v0*math.cos(math.radians(gamma0))
            p = h**2/planet.mu
            e = math.sqrt(1-p/a)
            ra = a*(1+e)
            hp = ((a*(1-e))- planet.Rp_e)/10**3  #(I'm considering a spherical planet)
            if hp < 0:
                print('WARNING AT initial_cond_calc: ALTITUDE PERIAPSIS < 0!')

            ra_list.append(ra)
            hp_list.append(hp)

    return ra_list, hp_list