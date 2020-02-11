from config import *
import math
import config

def nocontrol(m, rho, T_p, T_w, S, terminal_state):
    angle_of_attack = m.aerodynamics.aoa
    return angle_of_attack


def control_solarpanels_openloop(m,rho,T_p,T_w,S,terminal_state):
    taf  = m.aerodynamics.thermal_accomodation_factor
    R = m.planet.R
    gamma = m.planet.gamma
    max_angleofattack = math.pi/2
    min_angleofattack = 0.02
    thermal_limit = m.aerodynamics.thermal_limit - 0.00001
    noise = 0.0001
    x_0 = math.pi/3

    epsilon = 100
    heat_rate_max = heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, max_angleofattack)
    heat_rate_min = heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, min_angleofattack)
    counter = 0
    if (heat_rate_max >= thermal_limit) and (heat_rate_min <= thermal_limit):  # W/cm^2
        if abs(heat_rate_max - thermal_limit) < abs(heat_rate_min - thermal_limit):  #Newton method is unable to find a solution since there are multiple ones. We need to provide a good initial guess
            x_0 = math.pi/3
        elif abs(heat_rate_max - thermal_limit) > abs(heat_rate_min - thermal_limit):
            x_0 =  math.pi/6
        while (epsilon > 0.0001) and (counter < 500): # if get stuck in a local minimum
            counter += 1
            L = (taf * rho * R* T_p) * ((R * T_p / (2 * math.pi)) ** 0.5) * 10 ** -4
            fx = L * ((S ** 2 + (gamma) / (gamma - 1) - (gamma + 1) / (2 * (gamma - 1)) * (T_w / T_p)) * (math.exp(-(S * math.sin(x_0)) ** 2) + (math.pi ** 0.5) * (S * math.sin(x_0)) * (
                        1 + math.erf(S * math.sin(x_0)))) - 1 / 2 * math.exp(-(S * math.sin(x_0)) ** 2)) - thermal_limit
            #fx_prime = L*((math.pi)**0.5*S*math.cos(x_0)*(math.erf(S*math.sin(x_0))+1)*(-(gamma+1)/(2*(gamma-1)) * T_w/T_p +gamma/(gamma-1) + S**2) + (S**2 * math.sin(x_0) * math.cos(x_0) * math.exp(- (S*math.sin(x_0))**2)))
            fx_prime = L * S * math.cos(x_0) * ((math.pi ** 0.5) * (S ** 2 + gamma / (gamma - 1) + (gamma + 1) / (2 * (gamma - 1)) * T_w / T_p) * (
                    1 + math.erf(S * math.sin(x_0))) + S * math.sin(x_0) * math.exp(-(S * math.sin(x_0)) ** 2))
            x_1 = x_0 - (fx / fx_prime)
            epsilon = abs(x_1 - x_0)
            angle_of_attack = x_1
            x_0 = x_1
            epsilon_prev = epsilon
            if counter==499:
                print('Check controller - first check')
    elif(heat_rate_max < thermal_limit):
        angle_of_attack = math.pi / 2
    elif(heat_rate_min > thermal_limit):
        angle_of_attack = 0.00001
    else:
        print('Check controller - second check')


    if (angle_of_attack > math.pi / 2) or (angle_of_attack < 0):
        print(math.degrees(angle_of_attack))
        angle_of_attack = 0
        #print('Warning: No solution for Newton Method for angle of attack - change minimum periapsis!')


    # Evaluation Energy - porta a zero il controllo se ra z # ESPEDIENTE THERMAL LOAD
    #if (config.count_numberofpassage != 1) and ((config.solution.orientation.oe[0][-1] * (1 + config.solution.orientation.oe[1][-1])) - terminal_state[1] <= 10000*10**3):
    #    angle_of_attack = 0
    return angle_of_attack

def heat_rate_calc(taf, rho, T_w, T_p, R, gamma, S, angle):
    return (taf * rho * R * T_p) * ((R * T_p / (2 * math.pi)) ** 0.5) * ((S ** 2 + (gamma) / (gamma - 1) - (gamma + 1) / (2 * (gamma - 1)) * (T_w / T_p)) * (math.exp(-(S * math.sin(angle)) ** 2) + (math.pi ** 0.5) * (S * math.sin(angle)) *(1 + math.erf(S * math.sin(angle)))) - 1 / 2 * math.exp(-(S * math.sin(angle)) ** 2)) * 10 ** -4