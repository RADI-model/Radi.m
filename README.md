# RADI v1

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.4680228-informational)](https://doi.org/10.5281/zenodo.4739205)

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
2) Indicate whether you start the simulation from zero, continue an existing simulation, or start from predefined initial conditions.
"rerun=0;" is for starting from zero. This is the default option. 
"rerun=1;" is for continuing an existing simulation, where variables are already present in MATLAB's workspace.
"rerun=2;" is for loading a set of savec initial conditions.
3) Run RADI see the time steps at which results are saved being displayed in the commande window. 

Example 1:
```matlab
IC_W2 % set environmental conditions to those corresponding to the Equatorial Pacific station "W2"
rerun = 0; % for a fresh start
RADI % run the model
```
You can stop the run at anytime (ctrl+C), to change an parameter, plot the results, etc. Then, to continue the simulation from where it was stopped:
```matlab
rerun = 1; % for continuing an existing simulation 
RADI % run the model
```

Example 2:
```matlab
IC_W2 % set environmental conditions to those corresponding to the Equatorial Pacific station "W2"
load('steady_state_runs\ic_W2.mat') % load the initial conditions from a steady-state run at the W2 station
rerun = 2; % switch to loading initial conditions mode
RADI % run the model
```

## Changing the times at which RADI saves the results
The numbers that are displayed in the command represent the times after the start of the simulation (in year). During the simulation, RADI does not save the results at every time step, but only at predefined times. The default option is that RADI saves the results every year, indicated on line 80 of the IC_XXX files setting up environmental conditions:
```matlab
time_saved_resolution=1; %[a]
```
However, this is not adequate to visualize changes occuring at shorter timescales, e.g, tides, seasons, etc. 
In order to save the results every day or so, write:
```matlab
time_saved_resolution=1/320; %[a]
```
Note that it has to be a number with a limited number of decimals. In order to save the results every hour or so, write:
```matlab
time_saved_resolution=1/6400; %[a]
```
etc.

## Implementing time-dependent solid fluxes
In order to implement solid fluxes that are varying with time, go to the section entitled "Time varying Variable solid fluxes or DBL", line 361 and 365, in RADI.m.
Comment them out to have a Foc and Fcalcite that are varying with seasons. Another function can be implemented, i.e., to simulate a solid fluxes that is only decreasing through time. 

## Implementing time-dependent diffusive boundary layer thickness
In order to implement a diffusive boundary layer thickness that is varying with time, go to the section entitled "Time varying Variable solid fluxes or DBL", line 366, in RADI.m.
Comment it out to have a DBL that is varying with tides. 


