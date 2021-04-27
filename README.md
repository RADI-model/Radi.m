# RADI v1

1-D Reaction-Advection-Diffusion-Irrigation (RADI) Diagenetic Sediment Module.

Source code by [O. Sulpis](https://deep-c.xyz/), [M.P. Humphreys](https://mvdh.xyz), [M.M. Wilhelmus](https://wilhelmusmtz.github.io/) and D. Carroll. 

Requirements for RADI: [CO2SYS](https://github.com/jamesorr/CO2SYS-MATLAB/blob/master/src/CO2SYS.m).
Requirements for IC_W2, IC_SM7 and IC_H9: `gsw_rho` from the [Gibbs-SeaWater Oceanographic Toolbox](http://www.teos-10.org/software.htm).

Theory and applications are described in the accompanying model description paper: 
https://gmd.copernicus.org/articles/XXXXX


## Instructions

First, indicate the environmental conditions corresponding to the station you wish to simulate. 
Three stations are already included: 1) Southern Pacific (IC_SM7), 2) North Atlantic (IC_H9) and 3) Equatorial Pacific (IC_W2). 

If running for the first time:
```matlab
IC_W29 % set environmental conditions
rerun = 0; % for a fresh start; 1 for a re-run
RADI % run the model
```

```matlab
delta_phi = -phiBeta.*(phi0 - phiInf).*exp(-phiBeta*depths);
% delta_phi(1) = 0; % don't do this
delta_phiS = -delta_phi;
delta_tort2 = -2*delta_phi./phi; % not used in Julia
delta_D_bio = -2*depths.*D_bio/lambda_b^2; % not used in Julia
```
