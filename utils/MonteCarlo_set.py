import config
import numpy as np
import csv
import os

def MonteCarlo_setting():
    MC = {}
    MC['N passages'] , MC['Duration'] , MC['Median Heat'] , MC['Periapsis min'] , MC[
        'Periapsis max'] , count = [] , [] , [] , [] , [] , 0
    return MC, count

def MonteCarlo_setting_passage(mc_index,simulation):
    heat_passage = []
    print('--> MC number' , mc_index)
    # if commented: filename = current time (comment if you don't need)
    simulation['Filename'] = simulation['Filename'] + '_nMC={}'.format(mc_index + 1)

    config.index_MonteCarlo += 1
    # Initialization
    config.altitudeperiapsis , config.max_heatrate = [] , []
    return simulation

def MonteCarlo_append(MC, state , count):
    # Save results
    MC['N passages'].append(config.solution.orientation.numberofpassage[-1])
    MC['Duration'].append(config.solution.orientation.time[-1])
    MC['Median Heat'].append(np.median(config.max_heatrate))
    MC['Periapsis min'].append(min(config.altitudeperiapsis))
    MC['Periapsis max'].append(max(config.altitudeperiapsis))
    heat_rate_max = max(config.solution.performance.heat_rate)
    if heat_rate_max > state['Heat Rate']:
        count += 1
    print('--> Count =' , count)


def MonteCarlo_save(state, simulation,MC):
    folder_name = simulation['Filename'][0:simulation['Filename'].find('_nMC')]
    filename = '/Users/giusyfalcone/Aerobraking_SA_project_results/' + folder_name + '/MC_results_control={}'.format(simulation['Control Model'])+'_ra={}'.format(int(state['Apoapsis']/10**3))+'_rp={0:.1f}'.format(state['Periapsis'])+'_hl={}'.format(state['Heat Rate'])+'.csv'
    os.makedirs(os.path.dirname(filename) , exist_ok=True)
    with open(filename , "w") as f:
        writer = csv.writer(f , delimiter=',' , quotechar='"' , quoting=csv.QUOTE_MINIMAL)
        writer.writerow(range(simulation['Monte Carlo size']))
        writer.writerow(MC['N passages'])
        writer.writerow(MC['Duration'])
        writer.writerow(MC['Median Heat'])
        writer.writerow(MC['Periapsis min'])
        writer.writerow(MC['Periapsis max'])