import csv
from config import solution as s

def save_cvs(filename):
    with open(filename, "w") as f:
        #writer = csv.writer(f)
        writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        fieldnames = ['time','year','month','day','hour','min','sec','number of passage', 'pos_ii[0]', 'pos_ii[1]',
                      'pos_ii[2]', 'vel_ii[0]', 'vel_ii[1]', 'vel_ii[2]',
                      'pos_ii_mag', 'vel_ii_mag', 'pos_pp[0]', 'pos_pp[1]', 'pos_pp[2]', 'pos_pp_mag', 'vel_pp[0]',
                      'vel_pp[1]', 'vel_pp[2]', 'vel_pp_mag', 'a', 'e', 'i', 'OMEGA', 'omega', 'vi',
                      'lat', 'lon', 'alt', 'gamma_ii', 'gamma_pp', 'h_ii[0]', 'h_ii[1]', 'h_ii[2]', 'h_pp[0]',
                      'h_pp[1]', 'h_pp[2]', 'h_ii_mag', 'h_pp_mag', 'uD[0]', 'uD[1]', 'uD[2]', 'uE[0]', 'uE[1]',
                      'uE[2]', 'uN[0]', 'uN[1]', 'uN[2]', 'vN', 'vE', 'azi_pp', 'rho', 'T', 'p', 'wind[0]', 'wind[1]',
                      'wind[2]', 'cL', 'cD', 'aoa', 'mass', 'heat_rate', 'heat_load', 'T_r', 'q', 'gravity_ii[0]', 'gravity_ii[1]',
                      'gravity_ii[2]', 'drag_pp[0]', 'drag_pp[1]', 'drag_pp[2]', 'drag_ii[0]', 'drag_ii[1]',
                      'drag_ii[2]', 'lift_pp[0]', 'lift_pp[1]', 'lift_pp[2]', 'lift_ii[0]', 'lift_ii[1]', 'lift_ii[2]',
                      'force_ii[0]', 'force_ii[1]', 'force_ii[2]', 'energy']


        writer.writerow(fieldnames)
        for index in range(0 , len(s.orientation.time)):
            writer.writerow(
                [s.orientation.time[index] ,s.orientation.year[index], s.orientation.month[index], s.orientation.day[index],
                 s.orientation.hour[index], s.orientation.min[index], s.orientation.second[index],
                 s.orientation.numberofpassage[index], s.orientation.pos_ii[0][index] , s.orientation.pos_ii[1][index] ,
                 s.orientation.pos_ii[2][index] , s.orientation.vel_ii[0][index] , s.orientation.vel_ii[1][index] ,
                 s.orientation.vel_ii[2][index] , s.orientation.pos_ii_mag[index] , s.orientation.vel_ii_mag[index] ,
                 s.orientation.pos_pp[0][index] , s.orientation.pos_pp[1][index] , s.orientation.pos_pp[2][index] ,
                 s.orientation.pos_pp_mag[index] , s.orientation.vel_pp[0][index] , s.orientation.vel_pp[1][index] ,
                 s.orientation.vel_pp[2][index] , s.orientation.vel_pp_mag[index] , s.orientation.oe[0][index] ,
                 s.orientation.oe[1][index] ,s.orientation.oe[2][index] , s.orientation.oe[3][index] , s.orientation.oe[4][index] ,
                 s.orientation.oe[5][index] , s.orientation.lat[index],s.orientation.lon[index],s.orientation.alt[index],
                 s.orientation.gamma_ii[index], s.orientation.gamma_pp[index], s.orientation.h_ii[0][index],s.orientation.h_ii[1][index],s.orientation.h_ii[2][index],
                 s.orientation.h_pp[0][index],s.orientation.h_pp[1][index],s.orientation.h_pp[2][index],s.orientation.h_ii_mag[index],s.orientation.h_pp_mag[index],
                 s.orientation.uD[0][index],s.orientation.uD[1][index],s.orientation.uD[2][index],s.orientation.uE[0][index],s.orientation.uE[1][index],s.orientation.uE[2][index],
                 s.orientation.uN[0][index] , s.orientation.uN[1][index] , s.orientation.uN[2][index] ,  s.orientation.vN[index] , s.orientation.vE[index],s.orientation.azi_pp[index],
                 s.physical_properties.rho[index],s.physical_properties.T[index],s.physical_properties.p[index],s.physical_properties.wind[0][index],
                 s.physical_properties.wind[1][index] ,s.physical_properties.wind[2][index],s.physical_properties.cL[index],s.physical_properties.cD[index],s.physical_properties.aoa[index],
                 s.performance.mass[index] ,s.performance.heat_rate[index],s.performance.heat_load[index],s.performance.T_r[index],s.performance.q[index],
                 s.forces.gravity_ii[0][index], s.forces.gravity_ii[1][index], s.forces.gravity_ii[2][index],  s.forces.drag_pp[0][index], s.forces.drag_pp[1][index],
                 s.forces.drag_pp[2][index] , s.forces.drag_ii[0][index], s.forces.drag_ii[1][index], s.forces.drag_ii[2][index],
                 s.forces.lift_pp[0][index], s.forces.lift_pp[1][index], s.forces.lift_pp[2][index],s.forces.lift_ii[0][index],
                 s.forces.lift_ii[1][index] , s.forces.lift_ii[2][index] , s.forces.force_ii[0][index] ,s.forces.force_ii[1][index] ,s.forces.force_ii[2][index] ,
                 s.forces.energy[index],])
            #print(cnf.solution.orientation.time)
        #writer.writerow(variable)

        #writer = csv.DictWriter(f, fieldnames=fieldnames)
        #writer.writerows(a)
