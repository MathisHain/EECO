function dx = CO2atm_ode(t,x,wK,setSScsh,setCO2)

%=======================
%Parameters declariation
%=======================
% >> DATA FROM TOGGELWEIER AND SARMIENTO 1985 <<
T       = 24.5e6*1027*(60*60*24*365);              %Overturning,Kg/yr (m3/yr*1027Kg/m3) 24.5
M       = 38.1e6*1027*(60*60*24*365);              %High-latitude exchange, Kg/yr (m3/yr*1027Kg/m3); 38.1e6*1027 
V_tot   = 1.3e18*1027;              %Total ocean volume, m3
A_tot   = 1.3e18/3790;              %Ocean surface area, m2
mix_h   = 250 ;                     %Depth of mixed layer high lat ocean m
mix_l   = 100;                      %Depth of mixed layer low lat ocean m
area_h  = .15*A_tot;                %Area high lat ocean m2
area_l  = .85*A_tot;                %Area low lat ocean m2 
V_h     = area_h*mix_h*1027;         %Mass of high-latitude box (Area high is 15%) Kg (m3*1027Kg/m3)
V_l     = area_l*mix_l*1027;        %Mass of low-latitude box (Area low is 85%) Kg (m3*1027Kg/m3)
V_d     = (V_tot - V_h - V_l)*1027;  %Mass of deep box (calculated by difference) (m3*1027Kg/m3)
V_a     =1.773*10^20;                % Mass of atmospheric box [mol] 
%V_a     = 8.69e18*1.25;             %Mass TROPOSPHERE (calculated) (m3*1.25Kg/m3)

KE_l    = 0.99;
KE_h    = 0.2;

%global ALK_l ALK_h
temp_l  = 273.15+15; %Low latitude temperature [K] %%+15
temp_h  = 273.15+2; %High latitude temperature [K] %%2.0
salinity_l = 34.0; %low latitude salinity [g/Kg]
salinity_h = 31.0; %high latitude salinity [g/Kg]
%Even though the unit of measurement are different from the one used in the
%equations below (mmol/Kg vs mol/m3) this should not be a problem, as this
%parameters enters only in the function co2_flux which give us a flux which
%is moles/m2*s, which is compatible with our units.

Rcp = 116.0; % carbon to phosphorus ratio (to convert export form phosphorus to carbon unit)
Rnp = 16; % nitrogen to phosphorus ratio

eps = 5.5; %isotope effect for nitrate assimilation

Kwind=2/3600/24 ; % [m/s] piston velocity for air seas gas exchange///corresponding to 8.3 cm hr-1 (relatively close to 10cm hr-1 given in the Sarmiento&Gruber's book but depending on wind speed 

%[H, pH, pCO2, H2CO2, HCO3, CO3, Ozd, CSH] = carb_solver(x(22)+273.15,x(25),x(7), x(19), 4000) %deep ocean carbonate chemistry

%================
% DEFINE YOUR SYSTEM
%================

dx=zeros(25,1);
dx(1)=(-T/V_l)*x(1)+(T-KE_l*T)/V_l*x(3);
dx(2)=(T-KE_h*T)/V_h*x(1)+(-T-M)/V_h*x(2)+(M-KE_h*M)/V_h*x(3);
dx(3)=T/V_d*KE_h*x(1)+(T+M)/V_d*x(2)+(-M+KE_h*M+KE_l*T)/V_d*x(3);
dx(4)=0;
dx(5)=T/V_l*x(7)-T/V_l*x(5)-T/V_l*KE_l*Rcp*x(3)-(co2_flux(Kwind, temp_l, salinity_l,x(5), x(17), x(8),wK)*area_l)*(60*60*24*365)/V_l;
dx(6)=T/V_h*x(5)-T/V_h*x(6)+M/V_h*x(7)-M/V_h*x(6)-T/V_h*KE_h*Rcp*x(1)-M/V_h*KE_h*Rcp*x(3)-(co2_flux(Kwind, temp_h, salinity_h,x(6), x(18), x(8),wK)*area_h)*(60*60*24*365)/V_h;
dx(7)=T/V_d*x(6)-T/V_d*x(7)+M/V_d*x(6)-M/V_d*x(7)+T/V_d*KE_l*Rcp*x(3)+T/V_d*KE_h*Rcp*x(1)+M/V_d*KE_h*Rcp*x(3);
if isnan(setCO2)
	%disp('Keeping ddiagnostic CO2')
	dx(8)=((co2_flux(Kwind, temp_l, salinity_l,x(5), x(17), x(8),wK))*area_l)*(60*60*24*365)/V_a + ((co2_flux(Kwind, temp_h, salinity_h, x(6), x(18), x(8),wK))*area_h)*(60*60*24*365)/V_a;
 else
 	%disp('Holding CO2 levels constant at initialised initial condition')
	dx(8)=0;
end
	
dx(9)=T/V_l*x(11)-T/V_l*x(9)-T/V_l*KE_l*Rnp*x(3);
dx(10)=T/V_h*x(9)-T/V_h*x(10)+M/V_h*x(11)-M/V_h*x(10)-T/V_h*KE_h*Rnp*x(1)-M/V_h*KE_h*Rnp*x(3);
dx(11)=T/V_d*x(10)-T/V_d*x(11)+M/V_d*x(10)-M/V_d*x(11)+T/V_d*KE_l*Rnp*x(3)+T/V_d*KE_h*Rnp*x(1)+M/V_d*KE_h*Rnp*x(3);
dx(12)=0;
dx(13)=T/V_l*x(15)-T/V_l*x(13)-T/V_l*KE_l*Rnp*x(3)*x(15)/x(11)*(1-(1-KE_l)^(1-(eps)/1000))/KE_l;
dx(14)=T/V_h*x(13)-T/V_h*x(14)+M/V_h*x(15)-M/V_h*x(14)-T/V_h*KE_h*Rnp*x(1)*x(13)/x(9)*(1-(1-KE_h)^(1-(eps)/1000))/KE_h-M/V_h*KE_h*Rnp*x(3)*x(15)/x(11)*(1-(1-KE_h)^(1-(eps)/1000))/KE_h;
dx(15)=T/V_d*x(14)-T/V_d*x(15)+M/V_d*x(14)-M/V_d*x(15)+T/V_d*KE_l*Rnp*x(3)*x(15)/x(11)*(1-(1-KE_l)^(1-(eps)/1000))/KE_l+T/V_d*KE_h*Rnp*x(1)*x(13)/x(9)*(1-(1-KE_h)^(1-(eps)/1000))/KE_h+M/V_d*KE_h*Rnp*x(3)*x(15)/x(11)*(1-(1-KE_h)^(1-(eps)/1000))/KE_h;
dx(16)=0;

dx(17)=(T)/V_l*x(19)-(T)/V_l*x(17); % needs CaCO3 and OM-acidity export fluxes
dx(18)=(T)/V_h*x(17)+(M)/V_h*x(19)-(T+M)/V_h*x(18); % needs M-acidity export fluxe
dx(19)=(T+M)/V_d*x(18)-(T+M)/V_d*x(19); % needs dissolution/respiration

dx(20)=0; % set LL surf temp
dx(21)=0; % set HL surf temp
dx(22)=(T+M)/V_d*x(21)-(T+M)/V_d*x(22); % deep Sal from exchange with surface,needs checking

dx(23)=0; % set LL surf Sal
dx(24)=0; % set LL surf Sal
dx(25)=(T+M)/V_d*x(24)-(T+M)/V_d*x(25); % deep Sal from exchange with surface,needs checking

% CaCO3 compensation

if isnan(setSScsh)
	%disp('Keeping ddiagnostic CO2')
	dx(8)=((co2_flux(Kwind, temp_l, salinity_l,x(5), x(17), x(8),wK))*area_l)*(60*60*24*365)/V_a + ((co2_flux(Kwind, temp_h, salinity_h, x(6), x(18), x(8),wK))*area_h)*(60*60*24*365)/V_a;
 else
 	%disp('Holding CO2 levels constant at initialised initial condition')
	dx(8)=0;
end

end 



