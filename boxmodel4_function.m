function finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,wK, SSCSH, SSCO2) 
%%% to use the degree of nutrient concentration (KE_h) in the high-latitude box and impose a mean initial concentration for the whole ocean as for DIC and Alk
%%% what represent SSCSH?

% x vector:
% 1 = PO4_ll; 2 = PO4_hl; 3 = PO4_d; 4 = DIC_ll; 5 =  DIC_hl; 6 = DIC_D; 
% 7 = pCO2_a; 8 = Alk_ll; 9 = ALk_hl; 10 = Alk_d; 11 = T_ll; 12 = T_hl
% 13 = T_d; 14 = S_ll; 15 = S_hl; 16 = S_d

%%%%%%%     MODEL MAIN TEMPLATE    %%%%%%%
%==============================================
%INITIAL CONDITIONS and time span of simulation
%===============================================
    
PO4_ll_ini = 2.15*10^-6;%(mol/kg)
PO4_hl_ini = 2.15*10^-6;%(mol/kg)
PO4_d_ini  = 2.15*10^-6;%(mol/kg)
DIC_ll_ini = DICmean;%(mol/kg)1924*10^-6%%%%3000*10^-6
DIC_hl_ini = DICmean;%(mol/kg)2149*10^-6%%%%3351*10^-6
DIC_D_ini  = DICmean;%(mol/kg)2240*10^-6%%%%3493*10^-6

if isnan(SSCO2)
	disp('Using default initial CO2=280ppm')
	 pCO2_a_ini = 280*10^-6 ;%(atm))
 else
 	disp('Using set CO2 level')
 	 pCO2_a_ini = SSCO2*10^-6 ;%(atm))
 end

ALK_ll_ini = ALKmean;
ALK_hl_ini = ALKmean;
ALK_d_ini  = ALKmean;

T_ll_ini = 273.15+25;
T_hl_ini = 273.15+0;
T_d_ini = 273.15+5;

S_ll_ini = 35;
S_hl_ini = 34.7;
S_d_ini = 34.7;
  
x0 = [PO4_ll_ini, PO4_hl_ini, PO4_d_ini, DIC_ll_ini, DIC_hl_ini,DIC_D_ini,pCO2_a_ini, ALK_ll_ini, ALK_hl_ini, ALK_d_ini...
    T_ll_ini, T_hl_ini, T_d_ini, S_ll_ini,S_hl_ini,S_d_ini];

%=================
%SOLVING THE ODE
%=================
	tspan = (0:1:5000); %1000 years of simulation
    [t,x] = ode45(@(t,x)CO2atm_ode(t,x,KE_h,wK,SSCSH,SSCO2),tspan,x0,[]);
	finalstate = x;

end


