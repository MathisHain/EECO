function finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,wK, SSCSH, SSCO2,Tfeedback,init_dT,tmax) 
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

if isnan(SSCO2) & wK == "modern"
	disp('Using default initial CO2=280ppm')
	pCO2_a_ini = 280*10^-6 ;%(atm))
elseif isnan(SSCO2) & wK == "Eocene"
	disp('Using Eocene initial CO2=800ppm')
	pCO2_a_ini = 800*10^-6 ;%(atm))
 else
 	disp('Using set CO2 level')
 	pCO2_a_ini = SSCO2*10^-6 ;%(atm))
 end
 
 if isnan(SSCSH)
 	disp('No CaCO3 compensation; closed-system')
  else
  	fprintf('CaCO3 compensation with setCSH=%d; open-system\n',SSCSH)
  end

ALK_ll_ini = ALKmean;
ALK_hl_ini = ALKmean;
ALK_d_ini  = ALKmean;

T_ll_ini = 273.15+25+init_dT;
T_hl_ini = 273.15+0+init_dT;
T_d_ini = 273.15+5+init_dT;

S_ll_ini = 35;
S_hl_ini = 34.7;
S_d_ini = 34.7;
  
x0 = [PO4_ll_ini, PO4_hl_ini, PO4_d_ini, DIC_ll_ini, DIC_hl_ini,DIC_D_ini,pCO2_a_ini, ALK_ll_ini, ALK_hl_ini, ALK_d_ini...
    T_ll_ini, T_hl_ini, T_d_ini, S_ll_ini,S_hl_ini,S_d_ini];

%=================
%SOLVING THE ODE
%=================
	tspan = (0:1:tmax); %1000 years of simulation
    [t,x] = ode45(@(t,x)CO2atm_ode(t,x,KE_h,wK,SSCSH,SSCO2,Tfeedback),tspan,x0,[]);
	finalstate = x(1:100:end,:);

	V_tot   = 1.3e18*1027;              %Total ocean volume, m3
	A_tot   = 1.3e18/3790;              %Ocean surface area, m2
	mix_h   = 250 ;                     %Depth of mixed layer high lat ocean m
	mix_l   = 100;                      %Depth of mixed layer low lat ocean m
	area_h  = .15*A_tot;                %Area high lat ocean m2
	area_l  = .85*A_tot;                %Area low lat ocean m2 
	V_h     = area_h*mix_h*1027;         %Mass of high-latitude box (Area high is 15%) Kg (m3*1027Kg/m3)
	V_l     = area_l*mix_l*1027;        %Mass of low-latitude box (Area low is 85%) Kg (m3*1027Kg/m3)
	V_d     = (V_tot - V_h - V_l);
	ALKmean = (V_l*finalstate(:,8) + V_h*finalstate(:,9) + V_d*finalstate(:,10))/V_tot;
	DICmean = (V_l*finalstate(:,4) + V_h*finalstate(:,5) + V_d*finalstate(:,6))/V_tot;

	N = length(finalstate(:,8));
	CSH = zeros(N,1);
	for id = 1:N
		CSH(id) = carb_solver(finalstate(id,13),finalstate(id,16),finalstate(id,6), finalstate(id,10), 3000,wK);
	end	
	finalstate = [finalstate CSH ALKmean DICmean];
	
end


