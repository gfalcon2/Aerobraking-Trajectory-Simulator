import subprocess
# def call_ABTS(state):

from subprocess import Popen, PIPE

subprocess.call(["python3", "ABTS/ABTS.py", '--machine', 'Laptop','--integrator','Python', '--type_of_mission', 'Orbits','--number_of_orbits',str(1),
                '--control_mode',str(3),'--gravity_model','Inverse Squared',
                '--hp_initial_a', str(88000),'--density_model', 'MARSGram',
                 '--ra_initial_a', str(28038000),'--year', str(2001),'--month', str(12), '--day',str(14),'--hours', str(14),'--minutes', str(21), '--secs',str(28),
                 '--hp_step',str(10000000),'--ra_step',str(50000000000),'--MarsGram_version',str(0), '--max_heat_rate', str(0.15),'--max_heat_load', str(30),
               '--aop_dispersion_gnc', str(0), '--vi_dispersion_gnc', str(0),'--final_apoapsis',str(5088116.837416616 ),'--flash2_through_integration',str(0),'--flash1_rate',str(3),
                 '--control_in_loop',str(0),'--second_switch_reevaluation',str(1),'--security_mode',str(0)])

