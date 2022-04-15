#!/usr/bin/env python
from setuptools import setup, Extension
from Cython.Build import cythonize

import Cython.Compiler.Options
import numpy
Cython.Compiler.Options.annotate = True

ext_modules = cythonize([
    Extension("ABTS", ["ABTS.py"]) ,
    Extension("config", ["config.py"]),
    Extension("control.Control", ["control/Control.py"]),
    Extension("control.eoms" , ["control/eoms.py"]) ,
    Extension("control.Propulsive_maneuvers" , ["control/Propulsive_maneuvers.py"]) ,
    Extension("control.heatload_control.second_tsw_calcs" , ["control/heatload_control/second_tsw_calcs.py"]) ,
    Extension("control.heatload_control.security_mode" , ["control/heatload_control/security_mode.py"]) ,
    Extension("control.heatload_control.time_switch_calcs" , ["control/heatload_control/time_switch_calcs.py"]) ,
    Extension("control.heatload_control.utils_timeswitch" , ["control/heatload_control/utils_timeswitch.py"]) ,
    Extension("integrator.Integrators" , ["integrator/Integrators.py"]) ,
    Extension("integrator.Events" , ["integrator/Events.py"]) ,
    Extension("physical_models.Aerodynamic_models" , ["physical_models/Aerodynamic_models.py"]) ,
    Extension("physical_models.Density_models" , ["physical_models/Density_models.py"]) ,
    Extension("physical_models.Gravity_models" , ["physical_models/Gravity_models.py"]) ,
    Extension("physical_models.MARSGram" , ["physical_models/MarsGram.py"]) ,
    Extension("physical_models.Mission" , ["physical_models/Mission.py"]) ,
    Extension("physical_models.MonteCarlo_perturbations" , ["physical_models/MonteCarlo_perturbations.py"]) ,
    Extension("physical_models.Planet_data" , ["physical_models/Planet_data.py"]) ,
    Extension("physical_models.Propulsive_maneuvers" , ["physical_models/Propulsive_maneuvers.py"]),
    Extension("physical_models.Thermal_models" , ["physical_models/Thermal_models.py"]) ,
    Extension("simulation.run", ["simulation/run.pyx"]),
    Extension("simulation.Complete_passage" , ["simulation/Complete_passage.pyx"]) ,
    Extension("simulation.Aerobraking" , ["simulation/Aerobraking.py"]) ,
    Extension("simulation.Set_and_run" , ["simulation/Set_and_run.py"]) ,
    Extension("utils.closed_form_solution" , ["utils/closed_form_solution.py"]) ,
    Extension("utils.define_mission" , ["utils/define_mission.py"]) ,
    Extension("utils.initial_cond_calc" , ["utils/initial_cond_calc.py"]) ,
    Extension("utils.misc" , ["utils/misc.py"]) ,
    Extension("utils.MonteCarlo_set" , ["utils/MonteCarlo_set.py"]) ,
    Extension("utils.Odyssey_maneuver_plan" , ["utils/Odyssey_maneuver_plan.py"]) ,
    Extension("utils.plots" , ["utils/plots.py"]) ,
    Extension("utils.Ref_system_conf" , ["utils/Ref_system_conf.py"]) ,
    Extension("utils.Reference_system" , ["utils/Reference_system.py"]) ,
    Extension("utils.save_cvs" , ["utils/save_cvs.py"]) ,
    Extension("utils.save_results" , ["utils/save_results.py"]) ,
], annotate=True)


setup(
    ext_modules=ext_modules,
    include_dirs=[numpy.get_include()]
)

