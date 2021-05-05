%% 1-D Reaction-Advection-Diffusion-Irrigation (RADI) Diagenetic Sediment Module
%% Source code by O. Sulpis, M.P. Humphreys, M. Wilhelmus and D. Carroll
%% Uses: CO2SYS

disp("RADI is running the following setup/station:")
disp(Station)
if rerun==1
    disp("it is a rerun: stop the run at any time (CTRL+C) to visualize the evolution of the system and run again")
    disp("years after start of simulation:")
elseif rerun==0
    disp("it is not a rerun: stop the run at any time (CTRL+C) to visualize the evolution of the system")
    disp("then, set rerun to '1' to continue the analysis")
    disp("years after start of simulation")
else
    disp('initial conditions loaded')
    disp("years after start of simulation")
end

tStart = tic;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% carbonate system initialization this is used only at the first time step to initialize the calc_pco2 program
CO2SYS_data = CO2SYS(dtalkw*1e6/rho_sw,dtCO2w*1e6/rho_sw,1,2,S,T,T,P,P,dtSiw*1e6/rho_sw,dtPO4w*1e6/rho_sw,1,10,1);
k1(1,1:ndepths) = CO2SYS_data(1,67);    %carbonic acid first dissociation constant
k2(1,1:ndepths) = CO2SYS_data(1,68);     %carbonic acid second dissociation constant
k1p(1,1:ndepths) = CO2SYS_data(1,75);     %phosphate constant 1
k2p(1,1:ndepths) = CO2SYS_data(1,76);      %phosphate constant 2
k3p(1,1:ndepths) = CO2SYS_data(1,77);     %phosphate constant 3
kb(1,1:ndepths) = CO2SYS_data(1,72);      %boron constant 
kw(1,1:ndepths) = CO2SYS_data(1,71);       %water dissociation constants
ksi(1,1:ndepths) = CO2SYS_data(1,78);       %silica constants
bt(1,1:ndepths) = CO2SYS_data(1,79);      %[umol/kg] total boron 
omegaC = CO2SYS_data(1,30);          %calcite saturation state
omegaA = CO2SYS_data(1,31);          %aragonite saturation state
co3 = CO2SYS_data(1,22) .* 10^-6;          %[mol/kg] CO3 
hco3 = CO2SYS_data(1,21) .* 10^-6;          %[mol/kg] HCO3
ph = CO2SYS_data(1,37) ;         %pH on the total scale
Ca_ini = dCaw ./ rho_sw;         %[mol/kg] Ca concentration 
fg(1,1:ndepths)=dtalkw./rho_sw-hco3-2*co3;    %sum of all alkalinity species that are not carbon
kspc = (co3 .* Ca_ini) ./ omegaC;        %[mol2/kg2] calcite in situ solubility
kspa = (co3 .* Ca_ini) ./ omegaA;      %[mol2/kg2] aragonite in situ solubility
ff(1,1:ndepths) = 1;      %random parameter needed for calc_pco2
H(1,1:ndepths) = 10^-ph;       %[mol/kg] H concentration first guess
clear co3 hco3 ph omegaC omegaA Ca_ini CO2SYS_data   
sit = dtSiw./rho_sw;        %[mol/kg] convert silica concentration
bt = bt .* 10^-6;        %[mol/kg] convert boron concentration

%% temperature dependent "free solution" diffusion coefficients
D_dO2=0.031558+0.001428*T;       %[m2/a] oxygen diffusion coefficient 
D_dtalk=0.015179+0.000795*T;       %[m2/a] approximted to bicarbonate diffusion coefficient 
D_dtCO2=0.015179+0.000795*T;       %[m2/a] approximted to bicarbonate diffusion coefficient 
D_dtNO3=0.030863+0.001153*T;      %[m2/a] nitrate diffusion coefficient 
D_dtSO4=0.015779+0.000712*T;      %[m2/a] sulfate diffusion coefficient 
D_dtPO4=0.009783+0.000513*T;      %[m2/a] phosphate diffusion coefficient 
D_dtNH4=0.030926+0.001225*T;     %[m2/a] ammonium diffusion coefficient 
D_dtH2S=0.028938+0.001314*T;    %[m2/a] hydrogen sulfide diffusion coefficient 
D_dMn=0.009625+0.000481*T;    %[m2/a] manganese diffusion coefficient 
D_dFe=0.010761+0.000466*T;    %[m2/a] iron diffusion coefficient 
D_dCa=0.011771+0.000529*T;     %[m2/a] calcium diffusion coefficient 

%% bioturbation (for solids)
D_bio_0=1e-6*2.32*(Foc*1e2)^0.85;      %[m2/a] surf bioturb coeff, Archer et al (2002)
lambda_b = 0.08;      %[m] characteristic bioturbation depth
D_bio=D_bio_0*exp(-(depths./lambda_b).^2).*((dO2w/1e-3)/((dO2w/1e-3)+20)); %[m2/a] bioturb coeff, Archer et al (2002)

%% irrigation (for solutes)
alpha_0=11*(atan(((5*Foc*1e2-400))/400)/pi+0.5)-0.9...
    +20*((dO2w/1e-3)/((dO2w/1e-3)+10))*exp(-(dO2w/1e-3)/10)*Foc*1e2/(Foc*1e2+30);    %[/a] from Archer et al (2002)
lambda_i=0.05;    %[m] characteristic irrigation depth
alpha=alpha_0.*exp(-(depths/lambda_i).^2);    %[/a] Archer et al (2002) the depth of 5 cm was changed

%% depth-dependent porosity and diffusion coefficient loss
delta_phi = -phiBeta.*(phi0 - phiInf).*exp(-phiBeta*depths); % depth-dependent porosity loss
delta_phiS = -delta_phi;   % depth-dependent solid fraction gain
delta_tort2 = -2*delta_phi./phi;   % depth-dependent tortuosity gain [not used in Julia]
delta_D_bio = -2*depths.*D_bio/lambda_b^2; % [not used in Julia]

% biodiffusion depth-attenuation: see Boudreau (1996); Fiadeiro and Veronis (1977)
Peh=w.*z_res./(2*D_bio);      %one half the cell Peclet number (Eq. 97 in Boudreau 1996)
% when Peh<<1, biodiffusion dominates, when Peh>>1, advection dominates
sigma=1./tanh(Peh)-1./(Peh);  %Eq. 96 in Boudreau 1996

%% organic matter degradation parameters
 KdO2=0.003; %[mol/m3] Monod constant from Soetaert et al. 1996 (GCA)
 KindO2=0.01; %[mol/m3] Monod inhibition constant from Soetaert et al. 1996 (GCA)
 KdtNO3=0.03; %[mol/m3] Monod constant from Soetaert et al. 1996 (GCA)
 KindtNO3=0.005; %[mol/m3] Monod inhibition constant from Soetaert et al. 1996 (GCA)
 KpMnO2=42.4; %[mol/m3] Monod constant from Van Cappellen and Wang 1996
 KinpMnO2=KpMnO2; %[mol/m3] Monod inhibition constant from Van Cappellen and Wang 1996
 KpFeOH3=265; %[mol/m3] Monod constant from Van Cappellen and Wang 1996 
 KinpFeOH3=KpFeOH3; %[mol/m3] Monod inhibition constant from Van Cappellen and Wang 1996
 KdtSO4=1.6; %[mol/m3] Monod constant from Van Cappellen and Wang 1996
 KindtSO4=KdtSO4; %[mol/m3] Monod inhibition constant from Van Cappellen and Wang 1996

kslow_0=1.3e-4 * (Foc*1e2)^0.85;    %[/a] tuned parameter, function from Archer et al (2002)
kslow=kslow_0.*ones(1,size(depths,2));    %[/a] 

kfast_0=1.5e-1 * (Foc*1e2)^0.85;    %[/a] tuned parameter, function from Archer et al (2002)
kfast=kfast_0.*ones(1,size(depths,2));    %[/a] 

%% redox reaction first order constants
kMnox=1e6; %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)
kFeox=1e6; %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)
kNHox=1e4; %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)
kSox=3e5; %[mol/m3/a] rate constant for the deep-sea from Boudreau (1996)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% R.A.D.I. main loop %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rerun == 0 %concentrations set to zero for solids, bottom water values for not-too-sensitive solutes
    dO2(1,1:ndepths)=dO2w;        
    dtalk(1,1:ndepths)=dtalkw;            
    dtCO2(1,1:ndepths)=dtCO2w;            
    dtNO3(1,1:ndepths)=dtNO3w;      
    dtSO4(1,1:ndepths)=dtSO4w;      
    dtPO4(1,1:ndepths)=dtPO4w;      
    dtNH4(1,1:ndepths)=dtNH4w;      
    dtH2S(1,1:ndepths)=dtH2Sw;      
    dCa(1,1:ndepths)=dCaw;  
    dFe(1,1:ndepths)=dFew;
    dMn(1,1:ndepths)=dMnw;      
    proc(1,1:ndepths)=6e2;
    psoc(1,1:ndepths)=3e2;             
    pfoc(1,1:ndepths)=0;        
    pFeOH3(1,1:ndepths)=0;    
    pMnO2(1,1:ndepths)=0;    
    pcalcite(1,1:ndepths)=4.5e3;
    paragonite(1,1:ndepths)=0;
    pclay(1,1:ndepths)=1e3;
    
    % variable saving
    i=1;
    idx=1;
    plot_number=0:t_length/stoptime*time_saved_resolution:t_length;  %we will save the variables every "time_saved_resolution" year
    plot_number=0:t_length/stoptime*time_saved_resolution:t_length;  %we will save the variables every "time_saved_resolution" year
    plot_number(1)=1;
    
elseif rerun==1 %if it is a rerun, initial conditions are concentrations from last time step
    dO2=dO2f(:,idx-1)';            %[mol/m3]
    dtalk=dtalkf(:,idx-1)';            %[mol/m3]
    dtCO2=dtCO2f(:,idx-1)';      %[mol/m3]
    dtNO3=dtNO3f(:,idx-1)';            %[mol/m3]
    dtSO4=dtSO4f(:,idx-1)';            %[mol/m3]
    dtPO4=dtPO4f(:,idx-1)';            %[mol/m3]
    dtNH4=dtNH4f(:,idx-1)';            %[mol/m3]
    dtH2S=dtH2Sf(:,idx-1)';            %[mol/m3]
    dCa=dCaf(:,idx-1)';            %[mol/m3]
    dFe=dFef(:,idx-1)';            %[mol/m3]
    dMn=dMnf(:,idx-1)';            %[mol/m3]
    proc=procf(:,idx-1)';                %[mol/m3]
    psoc=psocf(:,idx-1)';                %[mol/m3]
    pfoc=pfocf(:,idx-1)';                %[mol/m3]
    pFeOH3=pFeOH3f(:,idx-1)';                %[mol/m3]
    pMnO2=pMnO2f(:,idx-1)';                %[mol/m3]
    pcalcite=pcalcitef(:,idx-1)';            %[mol/m3]
    paragonite=paragonitef(:,idx-1)';            %[mol/m3]
    pclay=pclayf(:,idx-1)';            %[mol/m3]
    
    plot_number=0:t_length/stoptime*time_saved_resolution:t_length;  %we will save the variables every "time_saved_resolution" year
    i=plot_number(idx-1);
    
else
    % initial condition for solutes: bottom-water value
    dO2=dO2ic';                %[mol/m3]
    dtalk=dtalkic';                %[mol/m3]
    dtCO2=dtCO2ic';                %[mol/m3]
    dtNO3=dtNO3ic';            %[mol/m3]
    dtSO4=dtSO4ic';            %[mol/m3]
    dtPO4=dtPO4ic';            %[mol/m3]
    dtNH4=dtNH4ic';            %[mol/m3]
    dtH2S=dtH2Sic';            %[mol/m3]
    dFe=dFeic';            %[mol/m3]
    dMn=dMnic';            %[mol/m3]
    dCa=dCaic';            %[mol/m3]
    proc=procic';
    psoc=psocic';
    pfoc=pfocic';
    pFeOH3=pFeOH3ic';
    pMnO2=pMnO2ic';
    pcalcite=pcalciteic';
    paragonite=paragoniteic';
    pclay=pclayic';
   end 
    % variable saving
    i=1;
    idx=1;
    plot_number=0:t_length/stoptime*time_saved_resolution:t_length;  %we will save the variables every "time_saved_resolution" year
    plot_number(1)=1;    


%% short-cut transport variables
APPW=w - delta_D_bio - delta_phiS.* D_bio./ phiS;
DFF=(tort2.* delta_phi./ phi - delta_tort2)./ (tort2.^2);
DBF=phiS.*D_bio+D_bio.*(-delta_phi);
TR=(2*z_res.* (tort.^2)./ dbl);

%% Prepare for timestep calculations // MPH [v20]

% Set indices for depth-varying reactions
j = 2:ndepths-1;
jp1 = j + 1;
jm1 = j - 1;
z_res2 = z_res.^2;

% Subsample variables that do not change from timestep to timestep
alpha_j = alpha(j);
z_res_j = z_res(j);
z_res2_j = z_res2(j);
DFF_j = DFF(j);
tort2_j = tort2(j);
u_j = u(j);
APPW_j = APPW(j);
sigma_j = sigma(j);
D_bio_j = D_bio(j);
phi_j = phi(j);
phiS_j = phiS(j);

%% Begin timestep loop

if rerun==0

dO2f = NaN(ndepths, stoptime+1);
dtalkf = NaN(ndepths, stoptime+1);
dtCO2f = NaN(ndepths, stoptime+1);
dtNO3f = NaN(ndepths, stoptime+1);
dtSO4f = NaN(ndepths, stoptime+1);
dtPO4f = NaN(ndepths, stoptime+1);
dtNH4f = NaN(ndepths, stoptime+1);
dtH2Sf = NaN(ndepths, stoptime+1);
dFef = NaN(ndepths, stoptime+1);
dMnf = NaN(ndepths, stoptime+1);
dCaf = NaN(ndepths, stoptime+1);
procf = NaN(ndepths, stoptime+1);
psocf = NaN(ndepths, stoptime+1);
pfocf = NaN(ndepths, stoptime+1);
pFeOH3f = NaN(ndepths, stoptime+1);
MnO2f = NaN(ndepths, stoptime+1);
pcalcitef = NaN(ndepths, stoptime+1);
paragonitef = NaN(ndepths, stoptime+1);
pclayf = NaN(ndepths, stoptime+1);
end

for i=i:t_length-1


    %%  Compute solute diffusive fluxes through the sediment-water interface 
    %    F_O2i=(D_dO2*phi(1)*(dO2f(1,:)-dO2w)./0.2e-3)';
    %    F_DICi=(D_dtCO2*phi(1)*(dtCO2f(1,:)-dtCO2w)./0.2e-3)';
    %    F_TAi=(D_dtalk*phi(1)*(dtalkf(1,:)-dtalkw)./0.2e-3)';
    
       
    %% Organic matter respiration pathways
    fdO2=dO2./(KdO2+dO2);                   %from the code of Couture et al. (EST 2010), following Boudreau (1996)
    fdtNO3=dtNO3./(KdtNO3+dtNO3).*(KindO2./(KindO2+dO2));
    fpMnO2=pMnO2./(pMnO2+KpMnO2).*(KindtNO3./(KindtNO3+dtNO3)).*(KindO2./(KindO2+dO2));
    fpFeOH3=pFeOH3./(pFeOH3+KpFeOH3).*(KinpMnO2./(pMnO2+KinpMnO2)).*(KindtNO3./(KindtNO3+dtNO3)).*(KindO2./(KindO2+dO2));
    fdtSO4=dtSO4./(dtSO4+KdtSO4).*(KinpFeOH3./(pFeOH3+KinpFeOH3)).*(KinpMnO2./(pMnO2+KinpMnO2)).*(KindtNO3./(KindtNO3+dtNO3)).*(KindO2./(KindO2+dO2));
    fdCH4=(KindtSO4./(dtSO4+KindtSO4)).*(KinpFeOH3./(pFeOH3+KinpFeOH3)).*(KinpMnO2./(pMnO2+KinpMnO2)).*(KindtNO3./(KindtNO3+dtNO3)).*(KindO2./(KindO2+dO2));
    fox=fdO2+fdtNO3+fpMnO2+fpFeOH3+fdtSO4+fdCH4;
    
    %% Redox reaction rates
    Rs_o2 = psoc.*kslow.*fdO2; %degradation of slow-decay organic matter by oxygen
    Rf_o2 = pfoc.*kfast.*fdO2; %degradation of fast-decay organic matter by oxygen
    Rs_no3 = psoc.*kslow.*fdtNO3; %...by nitrate
    Rf_no3 = pfoc.*kfast.*fdtNO3;
    Rs_mno2 = psoc.*kslow.*fpMnO2; %... by manganese oxide
    Rf_mno2 = pfoc.*kfast.*fpMnO2;
    Rs_feoh3 = psoc.*kslow.*fpFeOH3; %by iron oxide
    Rf_feoh3 = pfoc.*kfast.*fpFeOH3;
    Rs_so4 = psoc.*kslow.*fdtSO4; %... by sulfate
    Rf_so4 = pfoc.*kfast.*fdtSO4;
    Rs_ch4 = psoc.*kslow.*fdCH4; %... by itself
    Rf_ch4 = pfoc.*kfast.*fdCH4;
    Rs_tot = psoc.*kslow.*fox; % total degradation rate of slow-decay organic matter
    Rf_tot = pfoc.*kfast.*fox;  % total degradation rate of fast-decay organic matter
    RFeox = kFeox.*dFe.*dO2; % oxidation of dissolved iron
    RMnox = kMnox.*dMn.*dO2; % oxidation of dissolved manganese 
    RSox = kSox.*dtH2S.*dO2; % oxidation of hydrogen sulfide
    RNHox = kNHox.*dtNH4.*dO2; % oxidation of ammonia
    
    %% calc_pco2: time efficient carbonate system solver    
    DIC_molPerKg = dtCO2 ./ rho_sw;     %convert DIC to [mol/kg]
    TA_molPerKg = dtalk ./ rho_sw;         %convert TA to [mol/kg]
    PO4_molPerKg = dtPO4 ./ rho_sw;   %convert PO4 to [mol/kg]
    
    %======================================================================
    % Extract from calc_pCO2.m - runs ~3x faster here than in a separate
    % function // MPH [v20]
    %
    % Original comments:
    %
    %calc_pCO2()
    %Example FORTRAN subroutine: solve carbonate system for pC02
    %M. Follows, T. Ito, S. Dutkiewicz
    %D. Carroll, 2019
    hg = H;
    bohg = bt.*kb./(hg + kb);
    siooh3g = sit.*ksi./(ksi + hg);
    denom = hg.*hg.*hg + k1p.*hg.*hg + k1p.*k2p.*hg + k1p.*k2p.*k3p;
    h3po4g = PO4_molPerKg.*hg.*hg.*hg./denom;
    hpo4g = PO4_molPerKg.*k1p.*k2p.*hg./denom;
    po4g = PO4_molPerKg.*k1p.*k2p.*k3p./denom;

    %estimate carbonate alkalinity
    fg = (-bohg - (kw./hg)) + hg - hpo4g - 2*po4g + h3po4g - siooh3g;
    cag = TA_molPerKg + fg;

    %improved estimate of hydrogen ion conc
    gamm = DIC_molPerKg./cag;
    dummy = (1-gamm).*(1-gamm).*k1.*k1 - 4*k1.*k2.*(1 - 2*gamm);
    H = 0.5*((gamm-1).*k1 + sqrt(dummy)); %[mol/kg] H concentration
    % end calc_pCO2.m extract // MPH [v20]
    %======================================================================
    
    co3_molPerKg = (DIC_molPerKg .* k1 .* k2)./ (H.^2 + (k1 .* H) + (k1 .* k2)); %[mol/kg] CO3 concentration
    co3 = co3_molPerKg .* rho_sw; %[mol m^-3] CO3 concentraiton
    
    %% CaCO3 reactions
    OmegaC = dCa.*co3./ (kspc.* rho_sw.^2); %[no unit] calcite saturation state
    OmegaA = dCa.*co3./ (kspa.* rho_sw.^2); %[no unit] aragonite saturation state
    
   % Aragonite dissolution rate from Dong et al. (2019) EPSL
       for jj=1:ndepths
          if OmegaA(1,jj)<=1 && OmegaA(1,jj)>0.835 %this is the omega value for which both laws are equals
          Rd_aragonite(1,jj)=paragonite(1,jj).*0.0038*((1-OmegaA(1,jj)).^0.13);
          elseif OmegaA(1,jj)<=0.835 %this is the omega value for which both laws are equals
          Rd_aragonite(1,jj)=paragonite(1,jj).*0.0420*((1-OmegaA(1,jj)).^1.46);    
          else            
          Rd_aragonite(1,jj)=0;    
          end
       end
        
    %Calcite dissolution from Naviaux et al. (2019) Marine Chemistry
    for jj=1:ndepths
          if OmegaC(1,jj)<=1 && OmegaC(1,jj)>0.8275 %this is the omega value for which both laws are equals
          Rd_calcite(1,jj)=pcalcite(1,jj).*0.0063*((1-OmegaC(1,jj)).^0.11);
          elseif OmegaC(1,jj)<=0.8275 %this is the omega value for which both laws are equals
          Rd_calcite(1,jj)=pcalcite(1,jj).*20*((1-OmegaC(1,jj)).^4.7);    
          else            
          Rd_calcite(1,jj)=0;    
          end
    end
    
    %Calcite precipitation rate from Zuddas and Mucci, GCA (1998)
    %normalized to the same surface area than for dissolution (4m2/g)
    for jj=1:ndepths
         if OmegaC(1,jj)>1
         Rp_calcite(1,jj)=0.4075*((OmegaC(1,jj)-1).^1.76);  %%%CONSIDERABLE UNCERTAINTY IN THE RATE CONSTANT
         else 
         Rp_calcite(1,jj)=0;
         end
    end
        
        %% Time varying Variable solid fluxes or DBL
 
    %Foc=0.1956977592284939+sin(2*pi*i*interval/1)*0.1956977592284939*0.5; %[mol/m2/a] flux of total organic carbon to the bottom 
    %Froc=Foc*0.03; %[mol/m2/a] flux of refractory organic carbon to the bottom 
    %Fsoc=Foc*0.27; %[mol/m2/a] flux of slow-decay organic carbon to the bottom 
    %Ffoc=Foc*0.7; %[mol/m2/a] flux of fast-decay organic carbon to the bottom  
    %Fcalcite=0.22+sin(2*pi*i*interval/1)*0.22*0.5;
    %dbl=1e-3+0.5e-3*sin(2*pi*i*interval*365.25*24./6);
    
    %Foc_saved(1,idx)=Foc;
    %dbl_saved(1,idx)=dbl; 
    %TR=(2*z_res.* (tort.^2)./ dbl);

    %% Calculate all reactions (14 species, units: [mol/m3/a])
    TotR_dO2 = - phiS./phi.*(Rs_o2 + Rf_o2) - 0.25.*RFeox - 0.5.*RMnox - 2.*RSox - 2.*RNHox;
    TotR_dtalk = + phiS./ phi.*(Rs_o2.*(RN./RC - RP./RC) + Rf_o2.*(RN./RC - RP./RC) + Rs_no3.*(0.8+RN./RC - RP./RC)...
        + Rf_no3.*(0.8+RN./RC - RP./RC) + Rs_mno2.*(4+RN./RC - RP./RC) + Rf_mno2.*(4+RN./RC - RP./RC)...
        + Rs_feoh3.*(8+RN./RC - RP./RC) + Rf_feoh3.*(8+RN./RC - RP./RC) + Rs_so4.*(1+RN./RC - RP./RC)...
        + Rf_so4.*(1+RN./RC - RP./RC) + Rs_ch4.*(RN./RC - RP./RC) + Rf_ch4.*(RN./RC - RP./RC) + 2.* (Rd_calcite...
        + Rd_aragonite - Rp_calcite)) - 2.* RFeox - 2.* RMnox - 2.* RSox - 2.* RNHox;
    TotR_dtCO2 = phiS./phi.*(Rs_o2 + Rf_o2 + Rs_no3 + Rf_no3 + Rs_mno2 + Rf_mno2 + Rs_feoh3 + Rf_feoh3...
        + Rs_so4 + Rf_so4 + Rs_ch4.*0.5 + Rf_ch4.*0.5 + (Rd_calcite + Rd_aragonite - Rp_calcite));
    TotR_dtNO3 = - phiS./phi.*0.8.*(Rs_no3 + Rf_no3) + RNHox; 
    TotR_dtSO4 = - phiS./phi.*0.5.*(Rs_so4 + Rf_so4) + RSox; 
    TotR_dtPO4 = phiS./phi.*(RP./RC).*(Rs_tot + Rf_tot);
    TotR_dtNH4 = phiS./phi.*(RN./RC).*(Rs_tot + Rf_tot) - RNHox;
    TotR_dtH2S = phiS./phi.*0.5.*(Rs_so4 + Rf_so4) - RSox;    
    TotR_dFe = phiS./phi.*4.*(Rs_feoh3 + Rf_feoh3) - RFeox;
    TotR_dMn = phiS./phi.*2.*(Rs_mno2 + Rf_mno2) - RMnox;
    TotR_dCa = phiS./ phi.* (Rd_calcite + Rd_aragonite - Rp_calcite); 
    TotR_psoc = - Rs_tot;
    TotR_pfoc = - Rf_tot;
    TotR_pFeOH3 = - 4.*(Rs_feoh3 + Rf_feoh3) + phi./phiS.*RFeox;
    TotR_pMnO2 = - 2.*(Rs_mno2 + Rf_mno2) + phi./phiS.*RMnox;
    TotR_pcalcite = -Rd_calcite + Rp_calcite; 
    TotR_paragonite = -Rd_aragonite;
    
    %% top boundary condition: prescribed solid fluxes and diffusive boundary layer control on solutes
    dO2_1 = dO2(1) + interval * ( D_dO2 / tort2(1) * (2*dO2(2) - 2*dO2(1) + TR(1) * (dO2w - dO2(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dO2.*DFF(1)) * -1 * TR(1) * ( dO2w - dO2(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dO2w - dO2(1) ) ... %irrigation
        + TotR_dO2(1)); %reaction
    
     dtalk_1 = dtalk(1) + interval * ( D_dtalk / tort2(1) * (2*dtalk(2) - 2*dtalk(1) + TR(1) * (dtalkw - dtalk(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtalk.*DFF(1)) * -1 * TR(1) * ( dtalkw - dtalk(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtalkw - dtalk(1) ) ... %irrigation
        + TotR_dtalk(1)); %reaction
    
    dtCO2_1 = dtCO2(1) + interval * ( D_dtCO2 / tort2(1) * (2*dtCO2(2) - 2*dtCO2(1) + TR(1) * (dtCO2w - dtCO2(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtCO2.*DFF(1)) * -1 * TR(1) * ( dtCO2w - dtCO2(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtCO2w - dtCO2(1) ) ... %irrigation
        + TotR_dtCO2(1)); %reaction
    
    dtNO3_1 = dtNO3(1) + interval * ( D_dtNO3 / tort2(1) * (2*dtNO3(2) - 2*dtNO3(1) + TR(1) * (dtNO3w - dtNO3(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtNO3.*DFF(1)) * -1 * TR(1) * ( dtNO3w - dtNO3(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtNO3w - dtNO3(1) ) ... %irrigation
        + TotR_dtNO3(1)); %reaction
    
    dtSO4_1 = dtSO4(1) + interval * ( D_dtSO4 / tort2(1) * (2*dtSO4(2) - 2*dtSO4(1) + TR(1) * (dtSO4w - dtSO4(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtSO4.*DFF(1)) * -1 * TR(1) * ( dtSO4w - dtSO4(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtSO4w - dtSO4(1) ) ... %irrigation
        + TotR_dtSO4(1)); %reaction
    
    dtPO4_1 = dtPO4(1) + interval * ( D_dtPO4 / tort2(1) * (2*dtPO4(2) - 2*dtPO4(1) + TR(1) * (dtPO4w - dtPO4(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtPO4.*DFF(1)) * -1 * TR(1) * ( dtPO4w - dtPO4(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtPO4w - dtPO4(1) ) ... %irrigation
        + TotR_dtPO4(1)); %reaction
    
    dtNH4_1 = dtNH4(1) + interval * ( D_dtNH4 / tort2(1) * (2*dtNH4(2) - 2*dtNH4(1) + TR(1) * (dtNH4w - dtNH4(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtNH4.*DFF(1)) * -1 * TR(1) * ( dtNH4w - dtNH4(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtNH4w - dtNH4(1) ) ... %irrigation
        + TotR_dtNH4(1)); %reaction
    
    dtH2S_1 = dtH2S(1) + interval * ( D_dtH2S / tort2(1) * (2*dtH2S(2) - 2*dtH2S(1) + TR(1) * (dtH2Sw - dtH2S(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dtH2S.*DFF(1)) * -1 * TR(1) * ( dtH2Sw - dtH2S(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dtH2Sw - dtH2S(1) ) ... %irrigation
        + TotR_dtH2S(1)); %reaction
    
    dFe_1 = dFe(1) + interval * ( D_dFe / tort2(1) * (2*dFe(2) - 2*dFe(1) + TR(1) * (dFew - dFe(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dFe.*DFF(1)) * -1 * TR(1) * ( dFew - dFe(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dFew - dFe(1) ) ... %irrigation
        + TotR_dFe(1)); %reaction
    
    dMn_1 = dMn(1) + interval * ( D_dMn / tort2(1) * (2*dMn(2) - 2*dMn(1) + TR(1) * (dMnw - dMn(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dMn.*DFF(1)) * -1 * TR(1) * ( dMnw - dMn(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dMnw - dMn(1) ) ... %irrigation
        + TotR_dMn(1)); %reaction

    dCa_1 = dCa(1) + interval * ( D_dCa / tort2(1) * (2*dCa(2) - 2*dCa(1) + TR(1) * (dCaw - dCa(1))) / (z_res(1)^2) ... %diffusion
        - (u(1) - D_dCa.*DFF(1)) * -1 * TR(1) * ( dCaw - dCa(1)) / (2*z_res(1)) ... %advection
        + alpha(1) * ( dCaw - dCa(1) ) ... %irrigation
        + TotR_dCa(1)); %reaction

    proc_1 = proc(1) + interval * ( D_bio(1) * ( 2 * proc(2) - 2 * proc(1) +... %diffusion
        2 * z_res(1) * (Froc - phiS(1) * w(1) * proc(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*proc(2) + 2*sigma(1)*proc(1) - ... %advection
        (1 + sigma(1))*(proc(2)+2*z_res(1)/D_bio(1)*(Froc/phiS(1)-w(1)*proc(1))))/(2*z_res(1))); ... %advection
            
    psoc_1 = psoc(1) + interval * ( D_bio(1) * ( 2 * psoc(2) - 2 * psoc(1) +... %diffusion
        2 * z_res(1) * (Fsoc - phiS(1) * w(1) * psoc(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*psoc(2) + 2*sigma(1)*psoc(1) - ... %advection
        (1 + sigma(1))*(psoc(2)+2*z_res(1)/D_bio(1)*(Fsoc/phiS(1)-w(1)*psoc(1))))/(2*z_res(1)) ... %advection
        +TotR_psoc(1)); %reaction
    
    pfoc_1 = pfoc(1) + interval * ( D_bio(1) * ( 2 * pfoc(2) - 2 * pfoc(1) +... %diffusion
        2 * z_res(1) * (Ffoc - phiS(1) * w(1) * pfoc(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*pfoc(2) + 2*sigma(1)*pfoc(1) - ... %advection
        (1 + sigma(1))*(pfoc(2)+2*z_res(1)/D_bio(1)*(Ffoc/phiS(1)-w(1)*pfoc(1))))/(2*z_res(1)) ... %advection
        +TotR_pfoc(1)); %reaction
    
     pFeOH3_1 = pFeOH3(1) + interval * ( D_bio(1) * ( 2 * pFeOH3(2) - 2 * pFeOH3(1) +... %diffusion
        2 * z_res(1) * (FFeOH3 - phiS(1) * w(1) * pFeOH3(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*pFeOH3(2) + 2*sigma(1)*pFeOH3(1) - ... %advection
        (1 + sigma(1))*(pFeOH3(2)+2*z_res(1)/D_bio(1)*(FFeOH3/phiS(1)-w(1)*pFeOH3(1))))/(2*z_res(1)) ... %advection
        +TotR_pFeOH3(1)); %reaction

    pMnO2_1 = pMnO2(1) + interval * ( D_bio(1) * ( 2 * pMnO2(2) - 2 * pMnO2(1) +... %diffusion
        2 * z_res(1) * (FMnO2 - phiS(1) * w(1) * pMnO2(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*pMnO2(2) + 2*sigma(1)*pMnO2(1) - ... %advection
        (1 + sigma(1))*(pMnO2(2)+2*z_res(1)/D_bio(1)*(FMnO2/phiS(1)-w(1)*pMnO2(1))))/(2*z_res(1)) ... %advection
        +TotR_pMnO2(1)); %reaction
 
    pcalcite_1 = pcalcite(1) + interval * ( D_bio(1) * ( 2 * pcalcite(2) - 2 * pcalcite(1) +... %diffusion
        2 * z_res(1) * (Fcalcite - phiS(1) * w(1) * pcalcite(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*pcalcite(2) + 2*sigma(1)*pcalcite(1) - ... %advection
        (1 + sigma(1))*(pcalcite(2)+2*z_res(1)/D_bio(1)*(Fcalcite/phiS(1)-w(1)*pcalcite(1))))/(2*z_res(1)) ... %advection
        +TotR_pcalcite(1)); %reaction

    paragonite_1 = paragonite(1) + interval * ( D_bio(1) * ( 2 * paragonite(2) - 2 * paragonite(1) +... %diffusion
        2 * z_res(1) * (Faragonite - phiS(1) * w(1) * paragonite(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*paragonite(2) + 2*sigma(1)*paragonite(1) - ... %advection
        (1 + sigma(1))*(paragonite(2)+2*z_res(1)/D_bio(1)*(Faragonite/phiS(1)-w(1)*paragonite(1))))/(2*z_res(1)) ... %advection
        +TotR_paragonite(1)); %reaction

    pclay_1 = pclay(1) + interval * ( D_bio(1) * ( 2 * pclay(2) - 2 * pclay(1) +... %diffusion
        2 * z_res(1) * (Fclay - phiS(1) * w(1) * pclay(1)) / (D_bio(1) * phiS(1)) ) / (z_res(1).^2) ...  %diffusion
        - APPW(1) * ((1 - sigma(1))*pclay(2) + 2*sigma(1)*pclay(1) - ... %advection
        (1 + sigma(1))*(pclay(2)+2*z_res(1)/D_bio(1)*(Fclay/phiS(1)-w(1)*pclay(1))))/(2*z_res(1))); %advection
  
    %% bottom boundary condition: gradients disappear
    dO2_z = dO2(ndepths) + interval * (D_dO2 / tort2(ndepths) * 2 * ((dO2(ndepths-1) - dO2(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dO2w - dO2(ndepths)) ... %irrigation
        +TotR_dO2(ndepths));

    dtalk_z = dtalk(ndepths) + interval * (D_dtalk / tort2(ndepths) * 2 * ((dtalk(ndepths-1) - dtalk(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtalkw - dtalk(ndepths)) ... %irrigation
        +TotR_dtalk(ndepths));
    
    dtCO2_z = dtCO2(ndepths) + interval * (D_dtCO2 / tort2(ndepths) * 2 * ((dtCO2(ndepths-1) - dtCO2(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtCO2w - dtCO2(ndepths)) ... %irrigation
        +TotR_dtCO2(ndepths));
    
    dtNO3_z = dtNO3(ndepths) + interval * (D_dtNO3 / tort2(ndepths) * 2 * ((dtNO3(ndepths-1) - dtNO3(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtNO3w - dtNO3(ndepths)) ... %irrigation
        +TotR_dtNO3(ndepths));
    
    dtSO4_z = dtSO4(ndepths) + interval * (D_dtSO4 / tort2(ndepths) * 2 * ((dtSO4(ndepths-1) - dtSO4(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtSO4w - dtSO4(ndepths)) ... %irrigation
        +TotR_dtSO4(ndepths));
    
    dtPO4_z = dtPO4(ndepths) + interval * (D_dtPO4 / tort2(ndepths) * 2 * ((dtPO4(ndepths-1) - dtPO4(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtPO4w - dtPO4(ndepths)) ... %irrigation
        +TotR_dtPO4(ndepths));
    
    dtNH4_z = dtNH4(ndepths) + interval * (D_dtNH4 / tort2(ndepths) * 2 * ((dtNH4(ndepths-1) - dtNH4(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtNH4w - dtNH4(ndepths)) ... %irrigation
        +TotR_dtNH4(ndepths));
    
    dtH2S_z = dtH2S(ndepths) + interval * (D_dtH2S / tort2(ndepths) * 2 * ((dtH2S(ndepths-1) - dtH2S(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dtH2Sw - dtH2S(ndepths)) ... %irrigation
        +TotR_dtH2S(ndepths));
    
    dFe_z = dFe(ndepths) + interval * (D_dFe / tort2(ndepths) * 2 * ((dFe(ndepths-1) - dFe(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dFew - dFe(ndepths)) ... %irrigation
        +TotR_dFe(ndepths));
    
    dMn_z = dMn(ndepths) + interval * (D_dMn / tort2(ndepths) * 2 * ((dMn(ndepths-1) - dMn(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dMnw - dMn(ndepths)) ... %irrigation
        +TotR_dMn(ndepths));

    dCa_z = dCa(ndepths) + interval * (D_dCa / tort2(ndepths) * 2 * ((dCa(ndepths-1) - dCa(ndepths)) / z_res(ndepths).^2) ...  %diffusion
        + alpha(ndepths) * (dCaw - dCa(ndepths)) ... %irrigation
        +TotR_dCa(ndepths));
  
    proc_z = proc(ndepths) + interval * (D_bio(ndepths) * 2 * ( (proc(ndepths-1) - proc(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*proc(ndepths-1) + sigma(ndepths)*proc(ndepths))/z_res(ndepths)); %advection
    
    psoc_z = psoc(ndepths) + interval * (D_bio(ndepths) * 2 * ( (psoc(ndepths-1) - psoc(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*psoc(ndepths-1) + sigma(ndepths)*psoc(ndepths))/z_res(ndepths)... %advection
        +TotR_psoc(ndepths));
    
    pfoc_z = pfoc(ndepths) + interval * (D_bio(ndepths) * 2 * ( (pfoc(ndepths-1) - pfoc(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*pfoc(ndepths-1) + sigma(ndepths)*pfoc(ndepths))/z_res(ndepths)... %advection
        +TotR_pfoc(ndepths));
    
    pFeOH3_z = pFeOH3(ndepths) + interval * (D_bio(ndepths) * 2 * ( (pFeOH3(ndepths-1) - pFeOH3(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*pFeOH3(ndepths-1) + sigma(ndepths)*pFeOH3(ndepths))/z_res(ndepths)... %advection
        +TotR_pFeOH3(ndepths));
    
    pMnO2_z = pMnO2(ndepths) + interval * (D_bio(ndepths) * 2 * ( (pMnO2(ndepths-1) - pMnO2(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*pMnO2(ndepths-1) + sigma(ndepths)*pMnO2(ndepths))/z_res(ndepths)... %advection
        +TotR_pMnO2(ndepths));
    
    pcalcite_z = pcalcite(ndepths) + interval * (D_bio(ndepths) * 2 * ( (pcalcite(ndepths-1) - pcalcite(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*pcalcite(ndepths-1) + sigma(ndepths)*pcalcite(ndepths))/z_res(ndepths)... %advection
        +TotR_pcalcite(ndepths));

    paragonite_z = paragonite(ndepths) + interval * (D_bio(ndepths) * 2 * ( (paragonite(ndepths-1) - paragonite(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*paragonite(ndepths-1) + sigma(ndepths)*paragonite(ndepths))/z_res(ndepths)... %advection
        +TotR_paragonite(ndepths));
    
    pclay_z = pclay(ndepths) + interval * (D_bio(ndepths) * 2 * ( (pclay(ndepths-1) - pclay(ndepths)) / z_res(ndepths).^2)... %diffusion
        - APPW(ndepths) * (-sigma(ndepths)*pclay(ndepths-1) + sigma(ndepths)*pclay(ndepths))/z_res(ndepths)); %advection
           
    %% all other depths
    
    % Oxygen
    dO2_j = dO2(j);
    dO2_jp1 = dO2(jp1);
    dO2_jm1 = dO2(jm1);
    dO2(j) = dO2_j + interval*(TotR_dO2(j) - ...
        (u_j - D_dO2*DFF_j).*(dO2_jp1 - dO2_jm1)./(2*z_res_j) + ...
        (D_dO2./tort2_j).*((dO2_jp1 - 2*dO2_j + dO2_jm1)./z_res2_j) + ...
        alpha_j.*(dO2w - dO2_j));
    
    % Total alkalinity
    dtalk_j = dtalk(j);
    dtalk_jp1 = dtalk(jp1);
    dtalk_jm1 = dtalk(jm1);
    dtalk(j) = dtalk_j + interval*(TotR_dtalk(j) - ...
        (u_j - D_dtalk*DFF_j).*(dtalk_jp1 - dtalk_jm1)./(2*z_res_j) + ...
        (D_dtalk./tort2_j).*((dtalk_jp1 - 2*dtalk_j + dtalk_jm1)./z_res2_j) + ...
        alpha_j.*(dtalkw - dtalk_j));

    % Dissolved inorganic carbon
    dtCO2_j = dtCO2(j);
    dtCO2_jp1 = dtCO2(jp1);
    dtCO2_jm1 = dtCO2(jm1);
    dtCO2(j) = dtCO2_j + interval*(TotR_dtCO2(j) - ...
        (u_j - D_dtCO2*DFF_j).*(dtCO2_jp1 - dtCO2_jm1)./(2*z_res_j) + ...
        (D_dtCO2./tort2_j).*((dtCO2_jp1 - 2*dtCO2_j + dtCO2_jm1)./z_res2_j) + ...
        alpha_j.*(dtCO2w - dtCO2_j));

    % Nitrate
    dtNO3_j = dtNO3(j);
    dtNO3_jp1 = dtNO3(jp1);
    dtNO3_jm1 = dtNO3(jm1);
    dtNO3(j) = dtNO3_j + interval*(TotR_dtNO3(j) - ...
        (u_j - D_dtNO3*DFF_j).*(dtNO3_jp1 - dtNO3_jm1)./(2*z_res_j) + ...
        (D_dtNO3./tort2_j).*((dtNO3_jp1 - 2*dtNO3_j + dtNO3_jm1)./z_res2_j) + ...
        alpha_j.*(dtNO3w - dtNO3_j));
    
    % Sulfate
    dtSO4_j = dtSO4(j);
    dtSO4_jp1 = dtSO4(jp1);
    dtSO4_jm1 = dtSO4(jm1);
    dtSO4(j) = dtSO4_j + interval*(TotR_dtSO4(j) - ...
        (u_j - D_dtSO4*DFF_j).*(dtSO4_jp1 - dtSO4_jm1)./(2*z_res_j) + ...
        (D_dtSO4./tort2_j).*((dtSO4_jp1 - 2*dtSO4_j + dtSO4_jm1)./z_res2_j) + ...
        alpha_j.*(dtSO4w - dtSO4_j));
    
    % Phosphate
    dtPO4_j = dtPO4(j);
    dtPO4_jp1 = dtPO4(jp1);
    dtPO4_jm1 = dtPO4(jm1);
    dtPO4(j) = dtPO4_j + interval*(TotR_dtPO4(j) - ...
        (u_j - D_dtPO4*DFF_j).*(dtPO4_jp1 - dtPO4_jm1)./(2*z_res_j) + ...
        (D_dtPO4./tort2_j).*((dtPO4_jp1 - 2*dtPO4_j + dtPO4_jm1)./z_res2_j) + ...
        alpha_j.*(dtPO4w - dtPO4_j));
    
    % Ammonia
    dtNH4_j = dtNH4(j);
    dtNH4_jp1 = dtNH4(jp1);
    dtNH4_jm1 = dtNH4(jm1);
    dtNH4(j) = dtNH4_j + interval*(TotR_dtNH4(j) - ...
        (u_j - D_dtNH4*DFF_j).*(dtNH4_jp1 - dtNH4_jm1)./(2*z_res_j) + ...
        (D_dtNH4./tort2_j).*((dtNH4_jp1 - 2*dtNH4_j + dtNH4_jm1)./z_res2_j) + ...
        alpha_j.*(dtNH4w - dtNH4_j));
    
    % Hydrogen sulfide
    dtH2S_j = dtH2S(j);
    dtH2S_jp1 = dtH2S(jp1);
    dtH2S_jm1 = dtH2S(jm1);
    dtH2S(j) = dtH2S_j + interval*(TotR_dtH2S(j) - ...
        (u_j - D_dtH2S*DFF_j).*(dtH2S_jp1 - dtH2S_jm1)./(2*z_res_j) + ...
        (D_dtH2S./tort2_j).*((dtH2S_jp1 - 2*dtH2S_j + dtH2S_jm1)./z_res2_j) + ...
        alpha_j.*(dtH2Sw - dtH2S_j));
    
    % Dissolved iron
    dFe_j = dFe(j);
    dFe_jp1 = dFe(jp1);
    dFe_jm1 = dFe(jm1);
    dFe(j) = dFe_j + interval*(TotR_dFe(j) - ...
        (u_j - D_dFe*DFF_j).*(dFe_jp1 - dFe_jm1)./(2*z_res_j) + ...
        (D_dFe./tort2_j).*((dFe_jp1 - 2*dFe_j + dFe_jm1)./z_res2_j) + ...
        alpha_j.*(dFew - dFe_j));
    
    %Dissolved manganese
    dMn_j = dMn(j);
    dMn_jp1 = dMn(jp1);
    dMn_jm1 = dMn(jm1);
    dMn(j) = dMn_j + interval*(TotR_dMn(j) - ...
        (u_j - D_dMn*DFF_j).*(dMn_jp1 - dMn_jm1)./(2*z_res_j) + ...
        (D_dMn./tort2_j).*((dMn_jp1 - 2*dMn_j + dMn_jm1)./z_res2_j) + ...
        alpha_j.*(dMnw - dMn_j));
    
    % Dissolved calcium
    dCa_j = dCa(j);
    dCa_jp1 = dCa(jp1);
    dCa_jm1 = dCa(jm1);
    dCa(j) = dCa_j + interval*(TotR_dCa(j) - ...
        (u_j - D_dCa*DFF_j).*(dCa_jp1 - dCa_jm1)./(2*z_res_j) + ...
        (D_dCa./tort2_j).*((dCa_jp1 - 2*dCa_j + dCa_jm1)./z_res2_j) + ...
        alpha_j.*(dCaw - dCa_j));
    
   % Refractory organic carbon
    proc_j = proc(j);
    proc_jp1 = proc(jp1);
    proc_jm1 = proc(jm1);
    proc(j) = proc_j + interval*(... 
        - APPW_j.*(((1 - sigma_j).*proc_jp1 + ...
        2*sigma_j.*proc_j - ...
        (1 + sigma_j).*proc_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((proc_jp1 - 2*proc_j + ...
        proc_jm1)./z_res2_j));

    % Slow decay organic carbon
    psoc_j = psoc(j);
    psoc_jp1 = psoc(jp1);
    psoc_jm1 = psoc(jm1);
    psoc(j) = psoc_j + interval*(TotR_psoc(j) - ...
        APPW_j.*(((1 - sigma_j).*psoc_jp1 + ...
        2*sigma_j.*psoc_j - ...
        (1 + sigma_j).*psoc_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((psoc_jp1 - 2*psoc_j + ...
        psoc_jm1)./z_res2_j));
    
    % Fast decay organic carbon
    pfoc_j = pfoc(j);
    pfoc_jp1 = pfoc(jp1);
    pfoc_jm1 = pfoc(jm1);
    pfoc(j) = pfoc_j + interval*(TotR_pfoc(j) - ...
        APPW_j.*(((1 - sigma_j).*pfoc_jp1 + ...
        2*sigma_j.*pfoc_j - ...
        (1 + sigma_j).*pfoc_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((pfoc_jp1 - 2*pfoc_j + ...
        pfoc_jm1)./z_res2_j));
    
            % Iron oxide
    pFeOH3_j = pFeOH3(j);
    pFeOH3_jp1 = pFeOH3(jp1);
    pFeOH3_jm1 = pFeOH3(jm1);
    pFeOH3(j) = pFeOH3_j + interval*(TotR_pFeOH3(j) - ...
        APPW_j.*(((1 - sigma_j).*pFeOH3_jp1 + ...
        2*sigma_j.*pFeOH3_j - ...
        (1 + sigma_j).*pFeOH3_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((pFeOH3_jp1 - 2*pFeOH3_j + ...
        pFeOH3_jm1)./z_res2_j));    
    
            % Manganese oxide
    pMnO2_j = pMnO2(j);
    pMnO2_jp1 = pMnO2(jp1);
    pMnO2_jm1 = pMnO2(jm1);
    pMnO2(j) = pMnO2_j + interval*(TotR_pMnO2(j) - ...
        APPW_j.*(((1 - sigma_j).*pMnO2_jp1 + ...
        2*sigma_j.*pMnO2_j - ...
        (1 + sigma_j).*pMnO2_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((pMnO2_jp1 - 2*pMnO2_j + ...
        pMnO2_jm1)./z_res2_j));    
    
                % Calcite
    pcalcite_j = pcalcite(j);
    pcalcite_jp1 = pcalcite(jp1);
    pcalcite_jm1 = pcalcite(jm1);
    pcalcite(j) = pcalcite_j + interval*(TotR_pcalcite(j) - ...
        APPW_j.*(((1 - sigma_j).*pcalcite_jp1 + ...
        2*sigma_j.*pcalcite_j - ...
        (1 + sigma_j).*pcalcite_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((pcalcite_jp1 - 2*pcalcite_j + ...
        pcalcite_jm1)./z_res2_j));    
    
                % Aragonite
    paragonite_j = paragonite(j);
    paragonite_jp1 = paragonite(jp1);
    paragonite_jm1 = paragonite(jm1);
    paragonite(j) = paragonite_j + interval*(TotR_paragonite(j) - ...
        APPW_j.*(((1 - sigma_j).*paragonite_jp1 + ...
        2*sigma_j.*paragonite_j - ...
        (1 + sigma_j).*paragonite_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((paragonite_jp1 - 2*paragonite_j + ...
        paragonite_jm1)./z_res2_j));    
    
                % Clay
    pclay_j = pclay(j);
    pclay_jp1 = pclay(jp1);
    pclay_jm1 = pclay(jm1);
    pclay(j) = pclay_j + interval*( - ...
        APPW_j.*(((1 - sigma_j).*pclay_jp1 + ...
        2*sigma_j.*pclay_j - ...
        (1 + sigma_j).*pclay_jm1)./(2*z_res_j)) + ...
        D_bio_j.*((pclay_jp1 - 2*pclay_j + ...
        pclay_jm1)./z_res2_j));    
   
    %% Set top and bottom conditions in arrays
    dO2(1) = dO2_1;
    dtalk(1) = dtalk_1; 
    dtCO2(1) = dtCO2_1; 
    dtNO3(1) = dtNO3_1;
    dtSO4(1) = dtSO4_1;
    dtPO4(1) = dtPO4_1;
    dtNH4(1) = dtNH4_1;
    dtH2S(1) = dtH2S_1;
    dFe(1) = dFe_1; 
    dMn(1) = dMn_1;
    dCa(1) = dCa_1;
    proc(1) = proc_1;
    psoc(1) = psoc_1;
    pfoc(1) = pfoc_1; 
    pFeOH3(1) = pFeOH3_1;
    pMnO2(1) = pMnO2_1; 
    pcalcite(1) = pcalcite_1;
    paragonite(1) = paragonite_1;
    pclay(1) = pclay_1;
    dO2(ndepths) = dO2_z;
    dtalk(ndepths) = dtalk_z; 
    dtCO2(ndepths) = dtCO2_z; 
    dtNO3(ndepths) = dtNO3_z;
    dtSO4(ndepths) = dtSO4_z;
    dtPO4(ndepths) = dtPO4_z;
    dtNH4(ndepths) = dtNH4_z;
    dtH2S(ndepths) = dtH2S_z;
    dFe(ndepths) = dFe_z; 
    dCa(ndepths) = dCa_z;
    dMn(ndepths) = dMn_z;
    proc(ndepths) = proc_z;
    psoc(ndepths) = psoc_z;
    pfoc(ndepths) = pfoc_z; 
    pFeOH3(ndepths) = pFeOH3_z;
    pMnO2(ndepths) = pMnO2_z; 
    pcalcite(ndepths) = pcalcite_z;
    paragonite(ndepths) = paragonite_z;
    pclay(ndepths) = pclay_z;
        
    %% set very small or negative concentration to zero
    dO2(dO2<0)=0;
    dtalk(dtalk<0)=0;
    dtCO2(dtCO2<0)=0;
    dtNO3(dtNO3<0)=0;
    dtSO4(dtSO4<0)=0; 
    dtPO4(dtPO4<0)=0;
    dtNH4(dtNH4<0)=0; 
    dtH2S(dtH2S<0)=0; 
    dFe(dFe<0)=0; 
    dMn(dMn<0)=0;    
    dCa(dCa<0)=0;
    proc(proc<0)=0;
    psoc(psoc<0)=0;
    pfoc(pfoc<0)=0;
    pFeOH3(pFeOH3<0)=0;
    pMnO2(pMnO2<0)=0;    
    pcalcite(pcalcite<0)=0;    
    paragonite(paragonite<0)=0;    
    pclay(pclay<0)=0;    
    
    %% save data every step
    %dO2f(:, i+1) = dO2;
    %dtalkf(:, i+1) = dtalk;
    %dtCO2f(:, i+1) = dtCO2; 
    %dtNO3f(:, i+1) = dtNO3;
    %dtSO4f(:, i+1) = dtSO4;
    %dtPO4f(:, i+1) = dtPO4;
    %dtNH4f(:, i+1) = dtNH4;
    %dtH2Sf(:, i+1) = dtH2S;
    %dFef(:, i+1) = dFe; 
    %dMnf(:, i+1) = dMn;
    %dCaf(:, i+1) = dCa;
    %procf(:, i+1) = proc;
    %psocf(:, i+1) = psoc;
    %pfocf(:, i+1) = pfoc; 
    %pFeOH3f(:, i+1) = pFeOH3;
    %pMnO2f(:, i+1) = pMnO2; 
    %pcalcitef(:, i+1) = pcalcite; 
    %paragonitef(:, i+1) = paragonite; 
    %pclayf(:, i+1) = pclay; 
    
     if i == plot_number(idx)
       disp(plot_number(idx)*interval)
       dO2f(:, idx) = dO2;
       dtalkf(:, idx) = dtalk;
       dtCO2f(:, idx) = dtCO2; 
       dtNO3f(:, idx) = dtNO3;
       dtSO4f(:, idx) = dtSO4;
       dtPO4f(:, idx) = dtPO4;
       dtNH4f(:, idx) = dtNH4;
       dtH2Sf(:, idx) = dtH2S;
       dFef(:, idx) = dFe; 
       dMnf(:, idx) = dMn;
       dCaf(:, idx) = dCa;
       procf(:, idx) = proc;
       psocf(:, idx) = psoc;
       pfocf(:, idx) = pfoc; 
       pFeOH3f(:, idx) = pFeOH3;
       pMnO2f(:, idx) = pMnO2; 
       pcalcitef(:, idx) = pcalcite; 
       paragonitef(:, idx) = paragonite; 
       pclayf(:, idx) = pclay;
       idx=idx+1;
     end  
end

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
