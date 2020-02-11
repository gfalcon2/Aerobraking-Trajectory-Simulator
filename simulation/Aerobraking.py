import config as cnf
import math
from simulation.Complete_passage import asim
from utils.Ref_system_conf import OE
import time
import numpy as np

def aerobraking(ip, m, terminal_state, simulation):

    #Initialization
    cnf.count_aerobraking = 0
    cnf.count_overcome_hr = 0
    cnf.save_index_heat = 0
    FinalState = True
    continue_campaign = True
    numberofpassage = 0
    time_0 = 0
    initial_state = m.initialcondition

    # Aerobraking Campaign
    #while numberofpassage<1:
    while continue_campaign and (FinalState):
        cnf.index_Mars_Gram_call = 0
        numberofpassage = 1 + numberofpassage
        print("--> START PASSAGE #{}".format(numberofpassage))
        t = time.time()

        # Run Simulation
        continue_campaign = asim(ip, m, time_0, terminal_state, initial_state, numberofpassage, simulation)

        # Orbital Elements Result
        a = cnf.solution.orientation.oe[0][-1]
        e = cnf.solution.orientation.oe[1][-1]
        i = cnf.solution.orientation.oe[2][-1]
        OMEGA = cnf.solution.orientation.oe[3][-1]
        omega = cnf.solution.orientation.oe[4][-1]
        vi = math.pi+0.01

        # New Initial State
        initial_state = OE(a, e, i, OMEGA, omega, vi)

        v = [cnf.solution.orientation.vel_ii[0][-1], cnf.solution.orientation.vel_ii[1][-1], cnf.solution.orientation.vel_ii[2][-1]]
        r = [cnf.solution.orientation.pos_ii[0][-1], cnf.solution.orientation.pos_ii[1][-1], cnf.solution.orientation.pos_ii[2][-1]]
        time_0 = cnf.solution.orientation.time[-1]
        r_a = a*(1+e)
        elapsed = time.time() - t


        print('Mars Gram called {} times'.format(cnf.index_Mars_Gram_call))
        print('Computational Time {}s'.format(elapsed))
        print("--> PASSAGE #{} COMPLETE".format(numberofpassage))
        #print('Energy', ((np.inner(v, v)) / 2 - m.planet.mu / (np.inner(r, r)) ** 0.5))

        if simulation['Campaign Size'] == numberofpassage:
            continue_campaign = False

        if r_a <= terminal_state[1]:
            FinalState = False
            print('Reached Final State! R_a = ',r_a,'km')
            print('Thermal Limit overcomed totally', cnf.count_overcome_hr, 'times!')
        print(' ')

