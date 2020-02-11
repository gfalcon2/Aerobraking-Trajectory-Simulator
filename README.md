# Aerobraking Trajectory Simulator (ABTS)
Aerobraking trajectory simulator for different planets (Earth, Mars, Venus), missions included Odyssey and the possibility to add in-plane trajectory control during the drag passage.

The simulator was implemented by the [Putnam Research Group](https://putnam.aerospace.illinois.edu/) at the University of Illinois at Urbana-Champaign. The simulator has been used for the publication [Aerobraking Trajectory Control Using Articulated Solar Panels](https://www.researchgate.net/publication/339181606_Preprint_AAS_19-682_AEROBRAKING_TRAJECTORY_CONTROL_USING_ARTICULATED_SOLAR_PANELS).

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
simulation['Thermal Model'] = 'Maxwellian Heat Transfer'
simulation['Control Model'] = 0 #0='No Control',1='Rotative Solar Panels'
simulation['Only DragPassage'] = True
 ### Density Model
 The density can be evaluated through constant law, exponential law and through the use of MarsGram2010. However, the latter required [MarsGram2010](https://software.nasa.gov/software/MFS-33158-1) executable. The simulator is set to run MarsGram online. Wind effect can be taken into account if MarsGram2010 is used.

### Aerodynamic Model and Thermal Model
The aerodynamic coefficients are Mach-Dependent and are evaluated through the [Flow of Rarefied Gases](https://books.google.com/books?hl=en&lr=&id=DIIrDgAAQBAJ&oi=fnd&pg=PP1&dq=rarefied+flow+schaaf+and+chambre&ots=PWLd04BJmj&sig=DaKV6gVakAuvKRgQDM3ZE9uFrdQ#v=onepage&q=rarefied%20flow%20schaaf%20and%20chambre&f=false) theory by Schaaf and Chambre. The thermal model follows the Maxwellian Heat Transfer theory from the above reference.
