function steadystate = boxmodel4_function(Npaz,ALKmean)

%%%%%%%     MODEL MAIN TEMPLATE    %%%%%%%
%==============================================
%INITIAL CONDITIONS and time span of simulation
%===============================================
    
PO4_ll_ini = 0.5*10^-6;%(mol/kg)
PO4_hl_ini = 1.48*10^-6;%(mol/kg)
PO4_d_ini  = 2.15*10^-6;%(mol/kg)
PO4_a_ini  = 0.0; 
DIC_ll_ini = 1924*10^-6;%(mol/kg)1924*10^-6%%%%3000*10^-6
DIC_hl_ini = 2149*10^-6;%(mol/kg)2149*10^-6%%%%3351*10^-6
DIC_D_ini  = 2240*10^-6;%(mol/kg)2240*10^-6%%%%3493*10^-6
pCO2_a_ini = 280*10^-6 ;%(atm)) 
NO3_ll_ini = 5*10^-6;%(mol/kg)
NO3_hl_ini = Npaz;%25*10^-6;%(mol/kg)
NO3_d_ini  = 34.7*10^-6;%(mol/kg)
NO3_a_ini  = 0.0;   
%global ALK_l ALK_h
%ALK_l   = 2322*10^-6; %Low latitude alkalinity [mol/kg] - 2322*10^-6
%ALK_h   = 2322*10^-6; %high latitude alkalinity [mol/kg] - 2322*10^-6

d15N_ll_ini        = 5.4; % d15N for delta values
d15N_hl_ini        = 5.4;
d15N_d_ini         = 5.4;
d15N_a_ini         = 5.4;
R15ref             = 0.0036782; % isotopic ratio of Air N2 (i.e., the international reference for N isotopes)

N15_ll_ini         =((((d15N_ll_ini)/1000)+1)*R15ref)*NO3_ll_ini; %15N for 15N concentration
N15_hl_ini         =((((d15N_hl_ini)/1000)+1)*R15ref)*NO3_hl_ini;
N15_d_ini          =((((d15N_d_ini)/1000)+1)*R15ref)*NO3_d_ini;
N15_a_ini          =((((d15N_a_ini)/1000)+1)*R15ref)*NO3_a_ini;

x0 = [PO4_ll_ini, PO4_hl_ini, PO4_d_ini, PO4_a_ini, DIC_ll_ini, DIC_hl_ini,DIC_D_ini,pCO2_a_ini,NO3_ll_ini,NO3_hl_ini,NO3_d_ini,NO3_a_ini,N15_ll_ini,N15_hl_ini,N15_d_ini,N15_a_ini];

% need to keep track of ALK also as a dynamic state ariable
ALK_ll_ini = ALKmean;
ALK_hl_ini = ALKmean;
ALK_d_ini  = ALKmean;
x0 = [x0 [ALK_ll_ini, ALK_hl_ini, ALK_d_ini]]; % appending ALK as separate state variables to initial state vector


%=================
%SOLVING THE ODE
%=================
	tspan = (0:1:5000); %1000 years of simulation
    [t,x] = ode45(@CO2atm_ode,tspan,x0,[]);
	steadystate = x(end,:)

end


