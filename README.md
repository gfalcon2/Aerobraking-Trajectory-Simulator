# Aerobraking Trajectory Simulator (ABTS)
The Aerobraking Trajectory Simulator is a Python-based modeling and simulation tool developed to assess flight performance for the aerobraking mission set. The tool is open-source for the benefit of the entire aerospace community. The simulator can be used to model an entire aerobraking campaign, a single orbit, or a single atmospheric passage around Earth, Mars, and Venus, while providing estimates of important flight performance parameters, including
dynamic pressure, heat rate, and energy depletion rate. The computational efficiency of the tool also enables uncertainty analysis through Monte Carlo techniques. The built-in modularity of the tool enables selection of physical models from a range of built-in options, including a method to efficiently incorporate complex atmospheric models, such as Mars-GRAM 2010.
Aerobraking trajectory simulator for different planets (Earth, Mars, Venus), missions included Odyssey and the possibility to add in-plane trajectory control during the drag passage.

The simulator was implemented by the [Putnam Research Group](https://putnam.aerospace.illinois.edu/) at the University of Illinois at Urbana-Champaign. The simulator has been used for the publication [Aerobraking Trajectory Control Using Articulated Solar Panels](https://www.researchgate.net/publication/339181606_Preprint_AAS_19-682_AEROBRAKING_TRAJECTORY_CONTROL_USING_ARTICULATED_SOLAR_PANELS), [Closed-Form Trajectory Solution for Shallow, High-Altitude Atmospheric Flight](https://www.researchgate.net/profile/Giusy-Falcone/publication/344347043_AAS_20-448_Closed-Form_Trajectory_Solution_For_Shallow_High-Altitude_Atmospheric_Flight/links/5f6ab2eda6fdcc0086346859/AAS-20-448-Closed-Form-Trajectory-Solution-For-Shallow-High-Altitude-Atmospheric-Flight.pdf), [Energy Depletion Guidance for Aerobraking Atmospheric Passes](https://arc.aiaa.org/doi/abs/10.2514/1.G006171), and [Deep Reinforcement Learning for Autonomous Aerobraking Maneuver Planning](https://arc.aiaa.org/doi/abs/10.2514/6.2022-2497). 

Also, an overall overview of the simulator has been published in this paper [Design and Development of an Aerobraking Trajectory Simulation Tool](https://arc.aiaa.org/doi/abs/10.2514/6.2021-1065).

The simulator can be used to simulate an entire aerobraking campaign, one entire orbit, or only a drag-passage around Earth, Mars and Venus and it provides results in terms of trajectory, performances, physical properties, and forces. The Only-Drag-Passage simulation starts and ends at 160 km of altitude.

## Initial Conditions
Initial and final conditions have to be defined in Main.py. The initial conditions are apoapsis radius, periapsis altitude, inclination, RAAN, AOP, TA, while, the final conditions are only in terms of final apoapsis. For this reason, the final conditions are achieved only if a complete aerobraking campaign is simulated. For a drag passage simulation, the initial conditions can also be set to velocity and flight-path angle. 

## Vehicle
Although ABTS simulates the trajectory of a material point, the vehicle parameters (Cd, Cl, alfa, mass...) need to be defined in setting_and_run.py. Currently, the vehicle has the characteristic of the Odyssey spacecraft.

## Physical Models
### Planet
The planet is considered a rotative oblate. Latitude, longitude, and altitude are defined accordingly.
 ### Gravity Model
 Gravity model can be set to follow the constant law, the inverse squared low and taking into account the J2 effect.

 ### Density Model
 The density can be evaluated through constant law, exponential law and through the use of MarsGram2010. However, the latter required [MarsGram2010](https://software.nasa.gov/software/MFS-33158-1) executable. The simulator is set to run MarsGram online. Wind effect can be taken into account if MarsGram2010 is used.

### Aerodynamic Model and Thermal Model
The aerodynamic coefficients are Mach-Dependent and are evaluated through the [Flow of Rarefied Gases](https://books.google.com/books?hl=en&lr=&id=DIIrDgAAQBAJ&oi=fnd&pg=PP1&dq=rarefied+flow+schaaf+and+chambre&ots=PWLd04BJmj&sig=DaKV6gVakAuvKRgQDM3ZE9uFrdQ#v=onepage&q=rarefied%20flow%20schaaf%20and%20chambre&f=false) theory by Schaaf and Chambre. The thermal model follows the Maxwellian Heat Transfer theory from the above reference.

## Control 
The simulator enables the use of two forms of control: in-space propulsive control and in-atmosphere drag modulation trajectory control. 

### Propulsive Control:
In the case of propulsive control, ABTS can simulate aerobraking maneuvers, which are propulsive maneuvers performed at apoapsis to lower or raise the periapsis location, or drag passage thrust maneuvers for deceleration purposes. Generally, ABTS could be used for any other kind of propulsive maneuvers after designing them. To this end, the equations of motion are already set up to evaluate the trajectory with or without thrust. Finally, the propulsive
maneuvers require the definition of the Δ𝑉, which represents the estimate of the total change in velocity before and after the maneuver.

### Drag-Modulation Trajectory Control:
The simulator also allows a drag-modulation trajectory control during atmospheric drag passages through the rotation of the solar panels, as described in [Energy Depletion Guidance for Aerobraking Atmospheric Passes](https://arc.aiaa.org/doi/abs/10.2514/1.G006171). This control can be enabled only if the vehicle is a spacecraft.

This trajectory control aims to deplete the maximum amount of energy when the spacecraft flies in the atmospheric portion of the orbit, without exceeding pre-set thermal loads. This control method can be run in three configurations. The first configuration (1) constrains the heat rate; the second configuration (2) constrains the heat load; the last and third configuration (3) constrains the heat rate and load.

# Compile
The simulator can be used without compiling it. However, to speed up the simulation, it is suggested to compile ABTS.
## Prerequisite:
Cython
## Steps:
1) Open the terminal and navigate into the ABTS directory. 
2) Run the command: python3 setup.py build_ext —-inplace

Note: Remember to change the directory in which you wish to save your results in the args 'directory_results' in ABTS.py before compiling.
