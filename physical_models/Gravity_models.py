from utils.Reference_system import *
import numpy as np
def gravity_const(pos_ii_mag, pos_ii,p):
    pos_ii_hat = pos_ii / pos_ii_mag
    gravity_ii_mag_const = p.g_ref
    g = gravity_ii_mag_const * (pos_ii_hat)
    return g

def gravity_invsquared (pos_ii_mag, pos_ii, vel_ii, p):
    mu = p.mu
    pos_ii_hat = pos_ii / pos_ii_mag
    gravity_ii_mag_spherical = -mu/ pos_ii_mag ** 2
    g = gravity_ii_mag_spherical * (pos_ii_hat)
    return g

def gravity_invsquaredandJ2effect (pos_ii_mag, pos_ii, vel_ii, p):
    mu = p.mu
    J2 = p.J2
    r_i = cartesian(pos_ii[0], pos_ii[1], pos_ii[2])
    v_i = cartesian(vel_ii[0], vel_ii[1], vel_ii[2])

    method = 1

    if method == 1: # http://control.asu.edu/Classes/MAE462/462Lecture13.pdf
        # calculation of J2 in RTN system and transformation from RTN to inertial
        [OE] = rvtoorbitalelement(r_i, v_i, p)
        pos_ii_hat = pos_ii / pos_ii_mag
        gravity_ii_mag_spherical = -mu/ pos_ii_mag ** 2
        u = OE.omega+OE.vi #latitude angle:
        J2_rtn = -3*p.mu*J2*p.Rp_e**2/pos_ii_mag**4* np.array([1/2-3/2*np.sin(OE.i)**2* np.sin(u)**2,
                                                           np.sin(OE.i) ** 2 * np.sin(u)*np.cos(u),
                                                            np.sin(OE.i)*np.cos(OE.i)*np.sin(u)])

        T_ijk = [[np.cos(OE.OMEGA) * np.cos(u) - np.sin(OE.OMEGA) * np.sin(u) * np.cos(OE.i),
                  np.sin(OE.OMEGA) * np.cos(u) + np.cos(OE.OMEGA) * np.sin(u) * np.cos(OE.i),
                  np.sin(u) * np.sin(OE.i)],
                 [-np.cos(OE.OMEGA) * np.sin(u) - np.sin(OE.OMEGA) * np.cos(u) * np.cos(OE.i),
                  -np.sin(OE.OMEGA) * np.sin(u) + np.cos(OE.OMEGA) * np.cos(u) * np.cos(OE.i),
                  np.cos(u) * np.sin(OE.i)],
                 [np.sin(OE.OMEGA) * np.sin(OE.i), -np.cos(OE.OMEGA) * np.sin(OE.i), np.cos(OE.i)]]
        J2_ii = np.inner(np.transpose(T_ijk),J2_rtn)
        g = gravity_ii_mag_spherical * (pos_ii_hat) + J2_ii

    elif method == 2: # https://link.springer.com/content/pdf/10.1007/978-981-10-2383-5_2.pdf
        # calculation of J2 in LVLH system and transformation from LVLH to inertial
        [OE] = rvtoorbitalelement(r_i, v_i, p)
        u = OE.omega+OE.vi        #latitude angle:
        T_ijk = [[np.cos(OE.OMEGA) * np.cos(u) - np.sin(OE.OMEGA) * np.sin(u) * np.cos(OE.i),
                  -np.cos(OE.OMEGA) * np.sin(u) - np.sin(OE.OMEGA) * np.cos(u) * np.cos(OE.i),
                  np.sin(OE.OMEGA) * np.sin(OE.i)],
                 [np.sin(OE.OMEGA) * np.cos(u) + np.cos(OE.OMEGA) * np.sin(u) * np.cos(OE.i),
                  -np.sin(OE.OMEGA) * np.sin(u) + np.cos(OE.OMEGA) * np.cos(u) * np.cos(OE.i),
                  -np.cos(OE.OMEGA) * np.sin(OE.i)],
                 [np.sin(u) * np.sin(OE.i), np.cos(u)*np.sin(OE.i),np.cos(OE.i)]]

        # transformation of r from ECI to LVLH
        r_LHLV = np.inner(np.transpose(T_ijk),pos_ii)
        r_LHLV_mag = np.linalg.norm(r_LHLV)
        kJ2 = 3/2 *J2*mu*p.Rp_e**2
        grad_U = [mu/r_LHLV_mag**2 + (kJ2/r_LHLV_mag**4)*(1-3*np.sin(OE.i)**2*np.sin(u)**2),
                  (kJ2/r_LHLV_mag**4)*np.sin(OE.i)**2*np.sin(2*u),
                  (kJ2 / r_LHLV_mag ** 4) * np.sin(2*OE.i) * np.sin(u)]
        g_LHLV = [-item for item in grad_U]
        # transformation of gravitational acceleration from LVLH to ECI
        g = np.inner(T_ijk,g_LHLV)

    elif method == 3: # https://www.vcalc.com/equation/?uuid=1e5aa6ea-95a3-11e7-9770-bc764e2038f2
        # Calculation of J2 in ECI system
        pos_ii_hat = pos_ii / pos_ii_mag
        gravity_ii_mag_spherical = -mu/ pos_ii_mag ** 2
        x = pos_ii[0]
        y = pos_ii[1]
        z = pos_ii[2]
        J2_ii =3/2*J2*mu*(p.Rp_e**2)/(pos_ii_mag ** 5)* np.array([(5*(z/pos_ii_mag)**2 -1)*x,
                                                        (5 * (z / pos_ii_mag) ** 2 - 1) * y,
                                                       (5 * (z / pos_ii_mag) ** 2 - 3) * z])
        g = gravity_ii_mag_spherical * (pos_ii_hat) + J2_ii

    return g

# I tried calculation of J2 in different ways, all of them provide always the same results. From one passage to the other, periapsis altitude doesn't change that much, because since omega and OMEGA are set to be 0, the periapsis always occurs to the equator.
# However, the J2 effects the omega and OMEGA rates.