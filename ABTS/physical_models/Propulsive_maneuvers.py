import math
import config
from scipy import optimize
def propulsion_ic_calcs(m, args, initial_state):
    delta_v = args.delta_v * (-math.cos(args.phi))
    if args.print_res:
        if round(math.degrees(args.phi)) == 0:
            print("LOWER MANEUVER!")
        elif round(math.degrees(args.phi)) == 180:
            print("RAISE MANEUVER!")
    v_exausted = m.engine.Isp*m.engine.g_e
    if len(config.solution.performance.mass) == 0:
        m_i = m.body.Mass
    else:
        m_i = config.solution.performance.mass[-1]
    m_f = m_i/math.exp(delta_v/v_exausted)
    delta_m = (m_i-m_f)
    T = m.engine.T
    delta_t = delta_m*v_exausted/T
    delta_t_half = delta_t/2
    vi = initial_state.vi
    e = initial_state.e
    a = initial_state.a

    # check initial condition are not in drag pass 300 km
    r = m.planet.Rp_e +300*10**3
    vi_300 = math.acos((a * (1 - e ** 2) - r) / (e*r))
    E_300 = 2 * math.atan(((1 - e) / (1 + e)) ** 0.5 * math.tan((vi_300) / 2))  # eccentric anomaly at altitude 300 km
    E_apo = 2 * math.atan(((1 - e) / (1 + e)) ** 0.5 * math.tan((math.pi) / 2))  # eccentric anomaly at apoapsis
    delta_t_300 = ((a ** 3 / m.planet.mu) ** 0.5 * ((E_apo - e * math.sin(E_apo)) - (E_300 - e * math.sin(E_300))))
    if delta_t_half > delta_t_300:
        delta_t_half = delta_t_300
        delta_t = delta_t_half*2
        delta_m = delta_t*T/v_exausted
        m_f = m_i-delta_m
        delta_v = v_exausted*math.log(m_f/m_i)
        args.delta_v = delta_v/ (-math.cos(args.phi))
        print('-- Thrust Maximum Time Exceeded - Thrust Time and Delta-v adjusted - NEW DELTA-V = {} --'.format(args.delta_v))
    if len(config.solution.performance.mass) == 0:
        # Inverse Problem of Kepler
        E_0 = 2 * math.atan(
            ((1 - e) / (1 + e)) ** 0.5 * math.tan((vi) / 2))  # eccentric anomaly
        M_0 = E_0 - e*math.sin(E_0)
        n = 1/((a**3/m.planet.mu)**0.5)
        M_e = n*delta_t_half + M_0

        if M_e < math.pi:
            x_0 = M_e+e/2
        else:
            x_0 = M_e-e/2

        def f ( E ):
            fx = E-e*math.sin(E) -M_e
            return fx  # W/cm^2
        vi_in = optimize.newton(f, x_0, fprime=lambda E: 1-e*math.cos(E))

        if abs(vi_in) >= 2*math.pi:
            temp = round(vi_in/(2*math.pi))
            vi_in -= temp*2*math.pi
        if vi_in < 0:
            vi_in += 2*math.pi
        if vi_in > math.pi:
            delta = vi_in-math.pi
            vi_in = math.pi - delta
        initial_state.vi = vi_in
        return initial_state
    else:
        t_fin = config.solution.orientation.time[-1]
        t_in_new = t_fin - delta_t_half

        dist_list = [abs(t_in_new-p) for p in config.solution.orientation.time]
        temp = sorted(set(dist_list))
        index_a = dist_list.index(temp[0])
        index_b = index_a - 1
        dist_a = dist_list[index_a]
        dist_b = dist_list[index_b]

        t_in_new = (config.solution.orientation.time[index_a]+(config.solution.orientation.time[index_b]-config.solution.orientation.time[index_a])/(dist_b+abs(dist_a)) * abs(dist_a))
        results = [0]*90
        results[0] = (config.solution.orientation.year[index_a]+(config.solution.orientation.year[index_b]-config.solution.orientation.year[index_a])/(dist_b+abs(dist_a)) * abs(dist_a))
        results[1] = (config.solution.orientation.month[index_a] + (
                config.solution.orientation.month[index_b] - config.solution.orientation.month[index_a]) / (
                              dist_b + abs(dist_a)) * abs(dist_a))
        results[2] = (config.solution.orientation.day[index_a] + (
                config.solution.orientation.day[index_b] - config.solution.orientation.day[index_a]) / (
                              dist_b + abs(dist_a)) * abs(dist_a))

        results[3] = (config.solution.orientation.hour[index_a] + (
                config.solution.orientation.hour[index_b] - config.solution.orientation.hour[index_a]) / (
                              dist_b + abs(dist_a)) * abs(dist_a))
        results[4] = (config.solution.orientation.min[index_a] + (
                config.solution.orientation.min[index_b] - config.solution.orientation.min[index_a]) / (
                              dist_b + abs(dist_a)) * abs(dist_a))
        results[5] = (config.solution.orientation.second[index_a] + (
                config.solution.orientation.second[index_b] - config.solution.orientation.second[index_a]) / (
                              dist_b + abs(dist_a)) * abs(dist_a))

        results[6] = (config.solution.orientation.numberofpassage[index_a] + (
                config.solution.orientation.numberofpassage[index_b] - config.solution.orientation.numberofpassage[index_a]) / (
                              dist_b + abs(dist_a)) * abs(dist_a))
        results[7] = (config.solution.orientation.pos_ii[0][index_a] + (
                config.solution.orientation.pos_ii[0][index_b] - config.solution.orientation.pos_ii[0][index_a]) / (
                              dist_b + abs(dist_a)) * abs(dist_a))
        results[8] = (config.solution.orientation.pos_ii[1][index_a] + (
                config.solution.orientation.pos_ii[1][index_b] - config.solution.orientation.pos_ii[1][index_a]) / (
                              dist_b + abs(dist_a)) * abs(dist_a))
        results[9] = (config.solution.orientation.pos_ii[2][index_a] + (
                config.solution.orientation.pos_ii[2][index_b] - config.solution.orientation.pos_ii[2][index_a]) / (
                              dist_b + abs(dist_a)) * abs(dist_a))

        results[10] = (config.solution.orientation.vel_ii[0][index_a] + (
                config.solution.orientation.vel_ii[0][index_b] - config.solution.orientation.vel_ii[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[11] = (config.solution.orientation.vel_ii[1][index_a] + (
                config.solution.orientation.vel_ii[1][index_b] - config.solution.orientation.vel_ii[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[12] = (config.solution.orientation.vel_ii[2][index_a] + (
                config.solution.orientation.vel_ii[2][index_b] - config.solution.orientation.vel_ii[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[13] = (config.solution.orientation.pos_ii_mag[index_a] + (
                config.solution.orientation.pos_ii_mag[index_b] - config.solution.orientation.pos_ii_mag[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[14] = (config.solution.orientation.vel_ii_mag[index_a] + (
                config.solution.orientation.vel_ii_mag[index_b] - config.solution.orientation.vel_ii_mag[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[15] = (config.solution.orientation.pos_pp[0][index_a] + (
                config.solution.orientation.pos_pp[0][index_b] - config.solution.orientation.pos_pp[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[16] = (config.solution.orientation.pos_pp[1][index_a] + (
                config.solution.orientation.pos_pp[1][index_b] - config.solution.orientation.pos_pp[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[17] = (config.solution.orientation.pos_pp[2][index_a] + (
                config.solution.orientation.pos_pp[2][index_b] - config.solution.orientation.pos_pp[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[18] = (config.solution.orientation.pos_pp_mag[index_a] + (
                config.solution.orientation.pos_pp_mag[index_b] - config.solution.orientation.pos_pp_mag[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[19] = (config.solution.orientation.vel_pp[0][index_a] + (
                config.solution.orientation.vel_pp[0][index_b] - config.solution.orientation.vel_pp[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[20] = (config.solution.orientation.vel_pp[1][index_a] + (
                config.solution.orientation.vel_pp[1][index_b] - config.solution.orientation.vel_pp[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[21] = (config.solution.orientation.vel_pp[2][index_a] + (
                config.solution.orientation.vel_pp[2][index_b] - config.solution.orientation.vel_pp[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[22] = (config.solution.orientation.vel_pp_mag[index_a] + (
                config.solution.orientation.vel_pp_mag[index_b] - config.solution.orientation.vel_pp_mag[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[23] = (config.solution.orientation.oe[0][index_a] + (
                config.solution.orientation.oe[0][index_b] - config.solution.orientation.oe[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[24] = (config.solution.orientation.oe[1][index_a] + (
                config.solution.orientation.oe[1][index_b] - config.solution.orientation.oe[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[25] = (config.solution.orientation.oe[2][index_a] + (
                config.solution.orientation.oe[2][index_b] - config.solution.orientation.oe[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[26] = (config.solution.orientation.oe[3][index_a] + (
                config.solution.orientation.oe[3][index_b] - config.solution.orientation.oe[3][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[27] = (config.solution.orientation.oe[4][index_a] + (
                config.solution.orientation.oe[4][index_b] - config.solution.orientation.oe[4][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[28] = (config.solution.orientation.oe[5][index_a] + (
                config.solution.orientation.oe[5][index_b] - config.solution.orientation.oe[5][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[29] = (config.solution.orientation.lat[index_a] + (
                config.solution.orientation.lat[index_b] - config.solution.orientation.lat[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[30] = (config.solution.orientation.lon[index_a] + (
                config.solution.orientation.lon[index_b] - config.solution.orientation.lon[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[31] = (config.solution.orientation.alt[index_a] + (
                config.solution.orientation.alt[index_b] - config.solution.orientation.alt[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[32] = (config.solution.orientation.gamma_ii[index_a] + (
                config.solution.orientation.gamma_ii[index_b] - config.solution.orientation.gamma_ii[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[33] = (config.solution.orientation.gamma_pp[index_a] + (
                config.solution.orientation.gamma_pp[index_b] - config.solution.orientation.gamma_pp[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[34] = (config.solution.orientation.h_ii[0][index_a] + (
                config.solution.orientation.h_ii[0][index_b] - config.solution.orientation.h_ii[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[35] = (config.solution.orientation.h_ii[1][index_a] + (
                config.solution.orientation.h_ii[1][index_b] - config.solution.orientation.h_ii[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[36] = (config.solution.orientation.h_ii[2][index_a] + (
                config.solution.orientation.h_ii[2][index_b] - config.solution.orientation.h_ii[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[37] = (config.solution.orientation.h_pp[0][index_a] + (
                config.solution.orientation.h_pp[0][index_b] - config.solution.orientation.h_pp[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[38] = (config.solution.orientation.h_pp[1][index_a] + (
                config.solution.orientation.h_pp[1][index_b] - config.solution.orientation.h_pp[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[39] = (config.solution.orientation.h_pp[2][index_a] + (
                config.solution.orientation.h_pp[2][index_b] - config.solution.orientation.h_pp[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[40] = (config.solution.orientation.h_ii_mag[index_a] + (
                config.solution.orientation.h_ii_mag[index_b] - config.solution.orientation.h_ii_mag[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[41] = (config.solution.orientation.h_pp_mag[index_a] + (
                config.solution.orientation.h_pp_mag[index_b] - config.solution.orientation.h_pp_mag[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[42] = (config.solution.orientation.uD[0][index_a] + (
                config.solution.orientation.uD[0][index_b] - config.solution.orientation.uD[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[43] = (config.solution.orientation.uD[1][index_a] + (
                config.solution.orientation.uD[1][index_b] - config.solution.orientation.uD[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[44] = (config.solution.orientation.uD[2][index_a] + (
                config.solution.orientation.uD[2][index_b] - config.solution.orientation.uD[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[45] = (config.solution.orientation.uE[0][index_a] + (
                config.solution.orientation.uE[0][index_b] - config.solution.orientation.uE[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[46] = (config.solution.orientation.uE[1][index_a] + (
                config.solution.orientation.uE[1][index_b] - config.solution.orientation.uE[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[47] = (config.solution.orientation.uE[2][index_a] + (
                config.solution.orientation.uE[2][index_b] - config.solution.orientation.uE[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[48] = (config.solution.orientation.uN[0][index_a] + (
                config.solution.orientation.uN[0][index_b] - config.solution.orientation.uN[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[49] = (config.solution.orientation.uN[1][index_a] + (
                config.solution.orientation.uN[1][index_b] - config.solution.orientation.uN[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[50] = (config.solution.orientation.uN[2][index_a] + (
                config.solution.orientation.uN[2][index_b] - config.solution.orientation.uN[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[51] = (config.solution.orientation.vN[index_a] + (
                config.solution.orientation.vN[index_b] - config.solution.orientation.vN[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[52] = (config.solution.orientation.vE[index_a] + (
                config.solution.orientation.vE[index_b] - config.solution.orientation.vE[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[53] = (config.solution.orientation.azi_pp[index_a] + (
                config.solution.orientation.azi_pp[index_b] - config.solution.orientation.azi_pp[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[54] = (config.solution.physical_properties.rho[index_a] + (
                config.solution.physical_properties.rho[index_b] - config.solution.physical_properties.rho[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[55] = (config.solution.physical_properties.T[index_a] + (
                config.solution.physical_properties.T[index_b] - config.solution.physical_properties.T[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[56] = (config.solution.physical_properties.p[index_a] + (
                config.solution.physical_properties.p[index_b] - config.solution.physical_properties.p[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[57] = (config.solution.physical_properties.wind[0][index_a] + (
                config.solution.physical_properties.wind[0][index_b] - config.solution.physical_properties.wind[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[58] = (config.solution.physical_properties.wind[1][index_a] + (
                config.solution.physical_properties.wind[1][index_b] - config.solution.physical_properties.wind[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[59] = (config.solution.physical_properties.wind[2][index_a] + (
                config.solution.physical_properties.wind[2][index_b] - config.solution.physical_properties.wind[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[60] = (config.solution.physical_properties.cL[index_a] + (
                config.solution.physical_properties.cL[index_b] - config.solution.physical_properties.cL[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[61] = (config.solution.physical_properties.cD[index_a] + (
                config.solution.physical_properties.cD[index_b] - config.solution.physical_properties.cD[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[62] = (config.solution.physical_properties.aoa[index_a] + (
                config.solution.physical_properties.aoa[index_b] - config.solution.physical_properties.aoa[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[63] = (config.solution.physical_properties.S[index_a] + (
                config.solution.physical_properties.S[index_b] - config.solution.physical_properties.S[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[64] = (config.solution.performance.mass[index_a] + (
                config.solution.performance.mass[index_b] - config.solution.performance.mass[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[65] = (config.solution.performance.heat_rate[index_a] + (
                config.solution.performance.heat_rate[index_b] - config.solution.performance.heat_rate[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[66] = (config.solution.performance.heat_load[index_a] + (
                config.solution.performance.heat_load[index_b] - config.solution.performance.heat_load[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[67] = (config.solution.performance.T_r[index_a] + (
                config.solution.performance.T_r[index_b] - config.solution.performance.T_r[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[68] = (config.solution.performance.q[index_a] + (
                config.solution.performance.q[index_b] - config.solution.performance.q[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[69] = (config.solution.forces.gravity_ii[0][index_a] + (
                config.solution.forces.gravity_ii[0][index_b] - config.solution.forces.gravity_ii[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[70] = (config.solution.forces.gravity_ii[1][index_a] + (
                config.solution.forces.gravity_ii[1][index_b] - config.solution.forces.gravity_ii[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[71] = (config.solution.forces.gravity_ii[2][index_a] + (
                config.solution.forces.gravity_ii[2][index_b] - config.solution.forces.gravity_ii[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[72] = (config.solution.forces.drag_pp[0][index_a] + (
                config.solution.forces.drag_pp[0][index_b] - config.solution.forces.drag_pp[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[73] = (config.solution.forces.drag_pp[1][index_a] + (
                config.solution.forces.drag_pp[1][index_b] - config.solution.forces.drag_pp[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[74] = (config.solution.forces.drag_pp[2][index_a] + (
                config.solution.forces.drag_pp[2][index_b] - config.solution.forces.drag_pp[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[75] = (config.solution.forces.drag_ii[0][index_a] + (
                config.solution.forces.drag_ii[0][index_b] - config.solution.forces.drag_ii[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[76] = (config.solution.forces.drag_ii[1][index_a] + (
                config.solution.forces.drag_ii[1][index_b] - config.solution.forces.drag_ii[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[77] = (config.solution.forces.drag_ii[2][index_a] + (
                config.solution.forces.drag_ii[2][index_b] - config.solution.forces.drag_ii[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[78] = (config.solution.forces.lift_pp[0][index_a] + (
                config.solution.forces.lift_pp[0][index_b] - config.solution.forces.lift_pp[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[79] = (config.solution.forces.lift_pp[1][index_a] + (
                config.solution.forces.lift_pp[1][index_b] - config.solution.forces.lift_pp[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[80] = (config.solution.forces.lift_pp[2][index_a] + (
                config.solution.forces.lift_pp[2][index_b] - config.solution.forces.lift_pp[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[81] = (config.solution.forces.lift_ii[0][index_a] + (
                config.solution.forces.lift_ii[0][index_b] - config.solution.forces.lift_ii[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[82] = (config.solution.forces.lift_ii[1][index_a] + (
                config.solution.forces.lift_ii[1][index_b] - config.solution.forces.lift_ii[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[83] = (config.solution.forces.lift_ii[2][index_a] + (
                config.solution.forces.lift_ii[2][index_b] - config.solution.forces.lift_ii[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[84] = (config.solution.forces.force_ii[0][index_a] + (
                config.solution.forces.force_ii[0][index_b] - config.solution.forces.force_ii[0][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[85] = (config.solution.forces.force_ii[1][index_a] + (
                config.solution.forces.force_ii[1][index_b] - config.solution.forces.force_ii[1][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[86] = (config.solution.forces.force_ii[2][index_a] + (
                config.solution.forces.force_ii[2][index_b] - config.solution.forces.force_ii[2][index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[87] = (config.solution.forces.energy[index_a] + (
                config.solution.forces.energy[index_b] - config.solution.forces.energy[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        results[88] = (config.solution.simulation.MC_seed[index_a] + (
                config.solution.simulation.MC_seed[index_b] - config.solution.simulation.MC_seed[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))
        results[89] = (config.solution.simulation.drag_passage[index_a] + (
                config.solution.simulation.drag_passage[index_b] - config.solution.simulation.drag_passage[index_a]) / (
                               dist_b + abs(dist_a)) * abs(dist_a))

        # Clear all the results after the index
        ## SAVE RESULTS GLOBAL
        # Orientation
        del config.solution.orientation.time[index_b+1:]
        del config.solution.orientation.year[index_b+1:]
        del config.solution.orientation.month[index_b+1:]
        del config.solution.orientation.day[index_b+1:]
        del config.solution.orientation.hour[index_b+1:]
        del config.solution.orientation.min[index_b+1:]
        del config.solution.orientation.second[index_b+1:]
        del config.solution.orientation.numberofpassage[index_b+1:]
        del config.solution.orientation.pos_ii[0][index_b+1:]
        del config.solution.orientation.pos_ii[1][index_b+1:]
        del config.solution.orientation.pos_ii[2][index_b+1:]
        del config.solution.orientation.vel_ii[0][index_b+1:]
        del config.solution.orientation.vel_ii[1][index_b+1:]
        del config.solution.orientation.vel_ii[2][index_b+1:]
        del config.solution.orientation.pos_ii_mag[index_b+1:]
        del config.solution.orientation.vel_ii_mag[index_b+1:]

        del config.solution.orientation.pos_pp[0][index_b+1:]
        del config.solution.orientation.pos_pp[1][index_b+1:]
        del config.solution.orientation.pos_pp[2][index_b+1:]
        del config.solution.orientation.pos_pp_mag[index_b+1:]
        del config.solution.orientation.vel_pp[0][index_b+1:]
        del config.solution.orientation.vel_pp[1][index_b+1:]
        del config.solution.orientation.vel_pp[2][index_b+1:]
        del config.solution.orientation.vel_pp_mag[index_b+1:]

        del config.solution.orientation.oe[0][index_b+1:]
        del config.solution.orientation.oe[1][index_b+1:]
        del config.solution.orientation.oe[2][index_b+1:]
        del config.solution.orientation.oe[3][index_b+1:]
        del config.solution.orientation.oe[4][index_b+1:]
        del config.solution.orientation.oe[5][index_b+1:]

        del config.solution.orientation.lat[index_b+1:]
        del config.solution.orientation.lon[index_b+1:]
        del config.solution.orientation.alt[index_b+1:]
        del config.solution.orientation.gamma_ii[index_b+1:]
        del config.solution.orientation.gamma_pp[index_b+1:]

        del config.solution.orientation.h_ii[0][index_b+1:]
        del config.solution.orientation.h_ii[1][index_b+1:]
        del config.solution.orientation.h_ii[2][index_b+1:]
        del config.solution.orientation.h_pp[0][index_b+1:]
        del config.solution.orientation.h_pp[1][index_b+1:]
        del config.solution.orientation.h_pp[2][index_b+1:]
        del config.solution.orientation.h_ii_mag[index_b+1:]
        del config.solution.orientation.h_pp_mag[index_b+1:]

        del config.solution.orientation.uD[0][index_b+1:]
        del config.solution.orientation.uD[1][index_b+1:]
        del config.solution.orientation.uD[2][index_b+1:]
        del config.solution.orientation.uE[0][index_b+1:]
        del config.solution.orientation.uE[1][index_b+1:]
        del config.solution.orientation.uE[2][index_b+1:]
        del config.solution.orientation.uN[0][index_b+1:]
        del config.solution.orientation.uN[1][index_b+1:]
        del config.solution.orientation.uN[2][index_b+1:]
        del config.solution.orientation.vN[index_b+1:]
        del config.solution.orientation.vE[index_b+1:]
        del config.solution.orientation.azi_pp[index_b+1:]

        # Physical Properties
        del config.solution.physical_properties.rho[index_b+1:]
        del config.solution.physical_properties.T[index_b+1:]
        del config.solution.physical_properties.p[index_b+1:]
        del config.solution.physical_properties.wind[0][index_b+1:]
        del config.solution.physical_properties.wind[1][index_b+1:]
        del config.solution.physical_properties.wind[2][index_b+1:]
        del config.solution.physical_properties.cL[index_b+1:]
        del config.solution.physical_properties.cD[index_b+1:]
        del config.solution.physical_properties.aoa[index_b+1:]
        del config.solution.physical_properties.S[index_b+1:]

        # Performances
        del config.solution.performance.mass[index_b+1:]
        del config.solution.performance.heat_rate[index_b+1:]
        del config.solution.performance.heat_load[index_b+1:]
        del config.solution.performance.T_r[index_b+1:]
        del config.solution.performance.q[index_b+1:]

        # Forces
        del config.solution.forces.gravity_ii[0][index_b+1:]
        del config.solution.forces.gravity_ii[1][index_b+1:]
        del config.solution.forces.gravity_ii[2][index_b+1:]
        del config.solution.forces.drag_pp[0][index_b+1:]
        del config.solution.forces.drag_pp[1][index_b+1:]
        del config.solution.forces.drag_pp[2][index_b+1:]
        del config.solution.forces.drag_ii[0][index_b+1:]
        del config.solution.forces.drag_ii[1][index_b+1:]
        del config.solution.forces.drag_ii[2][index_b+1:]
        del config.solution.forces.lift_pp[0][index_b+1:]
        del config.solution.forces.lift_pp[1][index_b+1:]
        del config.solution.forces.lift_pp[2][index_b+1:]
        del config.solution.forces.lift_ii[0][index_b+1:]
        del config.solution.forces.lift_ii[1][index_b+1:]
        del config.solution.forces.lift_ii[2][index_b+1:]
        del config.solution.forces.force_ii[0][index_b+1:]
        del config.solution.forces.force_ii[1][index_b+1:]
        del config.solution.forces.force_ii[2][index_b+1:]
        del config.solution.forces.energy[index_b+1:]

        # Simulation
        del config.solution.simulation.MC_seed[index_b+1:]
        del config.solution.simulation.drag_passage[index_b+1:]

        config.solution.orientation.time.append(t_in_new)
        config.solution.orientation.year.append(results[0])
        config.solution.orientation.month.append(results[1])
        config.solution.orientation.day.append(results[2])
        config.solution.orientation.hour.append(results[3])
        config.solution.orientation.min.append(results[4])
        config.solution.orientation.second.append(results[5])
        config.solution.orientation.numberofpassage.append(results[6])
        config.solution.orientation.pos_ii[0].append(results[7])
        config.solution.orientation.pos_ii[1].append(results[8])
        config.solution.orientation.pos_ii[2].append(results[9])
        config.solution.orientation.vel_ii[0].append(results[10])
        config.solution.orientation.vel_ii[1].append(results[11])
        config.solution.orientation.vel_ii[2].append(results[12])
        config.solution.orientation.pos_ii_mag.append(results[13])
        config.solution.orientation.vel_ii_mag.append(results[14])

        config.solution.orientation.pos_pp[0].append(results[15])
        config.solution.orientation.pos_pp[1].append(results[16])
        config.solution.orientation.pos_pp[2].append(results[17])
        config.solution.orientation.pos_pp_mag.append(results[18])
        config.solution.orientation.vel_pp[0].append(results[19])
        config.solution.orientation.vel_pp[1].append(results[20])
        config.solution.orientation.vel_pp[2].append(results[21])
        config.solution.orientation.vel_pp_mag.append(results[22])

        config.solution.orientation.oe[0].append(results[23])
        config.solution.orientation.oe[1].append(results[24])
        config.solution.orientation.oe[2].append(results[25])
        config.solution.orientation.oe[3].append(results[26])
        config.solution.orientation.oe[4].append(results[27])
        config.solution.orientation.oe[5].append(results[28])

        config.solution.orientation.lat.append(results[29])
        config.solution.orientation.lon.append(results[30])
        config.solution.orientation.alt.append(results[31])
        config.solution.orientation.gamma_ii.append(results[32])
        config.solution.orientation.gamma_pp.append(results[33])

        config.solution.orientation.h_ii[0].append(results[34])
        config.solution.orientation.h_ii[1].append(results[35])
        config.solution.orientation.h_ii[2].append(results[36])
        config.solution.orientation.h_pp[0].append(results[37])
        config.solution.orientation.h_pp[1].append(results[38])
        config.solution.orientation.h_pp[2].append(results[39])
        config.solution.orientation.h_ii_mag.append(results[40])
        config.solution.orientation.h_pp_mag.append(results[41])

        config.solution.orientation.uD[0].append(results[42])
        config.solution.orientation.uD[1].append(results[43])
        config.solution.orientation.uD[2].append(results[44])
        config.solution.orientation.uE[0].append(results[45])
        config.solution.orientation.uE[1].append(results[46])
        config.solution.orientation.uE[2].append(results[47])
        config.solution.orientation.uN[0].append(results[48])
        config.solution.orientation.uN[1].append(results[49])
        config.solution.orientation.uN[2].append(results[50])
        config.solution.orientation.vN.append(results[51])
        config.solution.orientation.vE.append(results[52])
        config.solution.orientation.azi_pp.append(results[53])

        # Physical Properties
        config.solution.physical_properties.rho.append(results[54])
        config.solution.physical_properties.T.append(results[55])
        config.solution.physical_properties.p.append(results[56])
        config.solution.physical_properties.wind[0].append(results[57])
        config.solution.physical_properties.wind[1].append(results[58])
        config.solution.physical_properties.wind[2].append(results[59])
        config.solution.physical_properties.cL.append(results[60])
        config.solution.physical_properties.cD.append(results[61])
        config.solution.physical_properties.aoa.append(results[62])
        config.solution.physical_properties.S.append(results[63])

        # Performances
        config.solution.performance.mass.append(results[64])
        config.solution.performance.heat_rate.append(results[65])
        config.solution.performance.heat_load.append(results[66])
        config.solution.performance.T_r.append(results[67])
        config.solution.performance.q.append(results[68])

        # Forces
        config.solution.forces.gravity_ii[0].append(results[69])
        config.solution.forces.gravity_ii[1].append(results[70])
        config.solution.forces.gravity_ii[2].append(results[71])
        config.solution.forces.drag_pp[0].append(results[72])
        config.solution.forces.drag_pp[1].append(results[73])
        config.solution.forces.drag_pp[2].append(results[74])
        config.solution.forces.drag_ii[0].append(results[75])
        config.solution.forces.drag_ii[1].append(results[76])
        config.solution.forces.drag_ii[2].append(results[77])
        config.solution.forces.lift_pp[0].append(results[78])
        config.solution.forces.lift_pp[1].append(results[79])
        config.solution.forces.lift_pp[2].append(results[80])
        config.solution.forces.lift_ii[0].append(results[81])
        config.solution.forces.lift_ii[1].append(results[82])
        config.solution.forces.lift_ii[2].append(results[83])
        config.solution.forces.force_ii[0].append(results[84])
        config.solution.forces.force_ii[1].append(results[85])
        config.solution.forces.force_ii[2].append(results[86])
        config.solution.forces.energy.append(results[87])

        # Simulation
        config.solution.simulation.MC_seed.append(results[88])
        config.solution.simulation.drag_passage.append(results[89])

        return initial_state


