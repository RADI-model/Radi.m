# RADI v1

1-D Reaction-Advection-Diffusion-Irrigation (RADI) Diagenetic Sediment Module.

Source code by [O. Sulpis](https://deep-c.xyz/), [M.P. Humphreys](https://mvdh.xyz), [M.M. Wilhelmus](https://wilhelmusmtz.github.io/) and D. Carroll. 

Requirements for RADI: [CO2SYS](https://github.com/jamesorr/CO2SYS-MATLAB/blob/master/src/CO2SYS.m).
Requirements for IC_W2, IC_SM7 and IC_H9: `gsw_rho` from the [Gibbs-SeaWater Oceanographic Toolbox](http://www.teos-10.org/software.htm).
(also provided here)

Theory and applications are described in the accompanying model description paper: 
https://gmd.copernicus.org/articles/XXXXX


## Instructions

1) Load the environmental conditions corresponding to the station you wish to simulate. 
Three stations are already included: 1) Southern Pacific (IC_SM7), 2) North Atlantic (IC_H9) and 3) Equatorial Pacific (IC_W2). 
2) Indicate whether you start the simulation from zero, continue an existing simulation, or start from predefined initial conditions
"rerun=0;" is for starting from zero
"rerun=1;" is for continuing an existing simulation, where variables are already present in MATLAB's workspace
"rerun=2;" is for loading a set of savec initial conditions.
3) Run RADI see the time steps at which results are saved being displayed in the commande window. 


```matlab
IC_W29 % set environmental conditions
rerun = 0; % for a fresh start; 1 for a re-run
RADI % run the model
```
