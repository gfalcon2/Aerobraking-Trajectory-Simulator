from simulation.Set_and_run import *
import time as timecode
from utils.initial_cond_calc import *
from utils.MonteCarlo_set import *


## Initialization
state = {}
simulation = {}
# Save Results
simulation['Results'] = True
simulation['Print'] = False
simulation['Plot'] = True

## Model Definition
simulation['Planet'] = 1 # Earth: 0, Mars:1, Venus:2
simulation['Gravity Model'] = 'Inverse Squared'#''Inverse Squared and J2 effect'#'Inverse squared'#'Inverse Squared and J2 effect'
simulation['Density Model'] = 'MARSGram'#''MARSGram' #Exponential #'Constant'
simulation['Wind'] = 1 # wind calculation only if density model is MARSGram
simulation['Aerodynamic Model'] = 'Mach-dependent'  # 'Cd and Cl Constant'
simulation['Thermal Model'] = 'Maxwellian Heat Transfer'
simulation['Control Model'] = 0 #0='No Control',1='Rotative Solar Panels'
simulation['Only DragPassage'] = True
simulation['IC v-r'] = False # What initial condition do you have? if apoapsis radius and periapsis altitude = False/ if flight-path angle and velocity = True
simulation['Monte Carlo'] = 1
simulation['Monte Carlo size'] = 2
simulation['Campaign Size'] = 1

# Time
state['Year'], state['Month'], state['Day'], state['Hour'], state['Min'], state['Second'] = 2001, 11, 6, 8, 30, 0

# Maximum Thermal Condition
state['Heat Rate'] = 0.1  #W/cm^2

# CASES:
# Case 1: only drag passage, no MC
# Case 2: only drag passage, yes MC
# Case 3: all campaign, no MC
# Case 4: all campaign, yes MC
# Case 5: User case

Case = 5
if Case == 1:
    simulation['Only DragPassage'] = True
    simulation['Monte Carlo'] = False
elif Case == 2:
    simulation['Only DragPassage'] = True
    simulation['Monte Carlo'] = True
elif Case == 3:
    simulation['Only DragPassage'] = False
    simulation['Monte Carlo'] = False
elif Case == 4:
    simulation['Only DragPassage'] = False
    simulation['Monte Carlo'] = True
else:
    pass



## SIMULATION
t = timecode.time()
if simulation['Only DragPassage'] == True:  ## The passage would be simulated only in the atmosphere, the campaign size level would be constricted at 1 and the results will be only about drag passage (h < 160 km)
    simulation['Campaign Size'] = 1

# Initial Condition # Define initial flight path angle and velocity
if (simulation['Only DragPassage'] == True) and (simulation['IC v-r'] == True):
    gamma_0 , v_0 = [-5 , -4] , [4500 , 4600]  # deg and m/s
    planet = planet_data(simulation)
    apoapsis , periapsis_alt = ic_calculation_rptoae(planet , gamma_0 , v_0)
    inclination , omega = 93.6 , 0

else:
    apoapsis , periapsis_alt , inclination , omega = [28000 * 10 ** 3] , [110.] ,93.6,0 #93.6 , 0 28523.95142557378  # with values of omega not equal to 0, recalculate periapsis altitude due to the oblateness of the planet
    final_apoapsis = 4905.974818462152 * 10 ** 3

    if simulation['Monte Carlo'] == False:
        simulation['Monte Carlo size'] = 1
    else:
        MC,count = MonteCarlo_setting()

    for periapsis_item in periapsis_alt:
        for apoapsis_item in apoapsis:

           # Orbital Element
            state['Apoapsis'] , state['Periapsis'] , state['Inclination'] , state['OMEGA'] , state[
                'omega'] , state[
                'Final SMA'] = apoapsis_item , periapsis_item , inclination , omega , omega , final_apoapsis


            for mc_index in range(simulation['Monte Carlo size']):
                simulation['Filename'] = 'Results_control={}'.format(simulation['Control Model']) + '_ra={}'.format(int(state['Apoapsis'] / 10 ** 3)) + '_rp={0:.1f}'.format(state['Periapsis']) + '_hl={}'.format(state['Heat Rate'])

                if simulation['Monte Carlo'] == True:
                    MonteCarlo_setting_passage(mc_index, simulation)

                aeroraking_campaign(state, simulation)

                if simulation['Monte Carlo'] == True:
                    MonteCarlo_append(MC , state , count)

            if simulation['Monte Carlo'] == True:
                MonteCarlo_save(state, simulation , MC)
    print('--> COMPUTATIONAL TIME = ' , timecode.time() - t)



