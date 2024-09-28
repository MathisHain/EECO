function dx = CO2atm_ode(t,x,KE_h,wK,setSScsh,setCO2)

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
V_d     = (V_tot - V_h - V_l);  	%Mass of deep box (calculated by difference) (m3*1027Kg/m3)
V_a     =1.773*10^20;                % Mass of atmospheric box [mol] 
%V_a     = 8.69e18*1.25;             %Mass TROPOSPHERE (calculated) (m3*1.25Kg/m3)

KE_l    = 0.99;
%KE_h    = 0.2;

%global ALK_l ALK_h
temp_l  = 273.15+x(11); %Low latitude temperature [K] 
temp_h  = 273.15+x(12); %High latitude temperature [K] 
salinity_l = x(14);%low latitude salinity [g/Kg]
salinity_h = x(15);%high latitude salinity [g/Kg]

Rpump = 0.2; % carbonate vs soft-tissue pump in low-latitude area
Rcp = 116.0; % carbon to phosphorus ratio (to convert export form phosphorus to carbon unit)
Rnp = 16; % nitrogen to phosphorus ratio (for the effect of the soft-tissue pump on alkalinity
Ralk = 2; % Alk to DIC ratio for the carbonate pump

Kwind=2/3600/24 ; % [m/s] piston velocity for air seas gas exchange///corresponding to 8.3 cm hr-1 (relatively close to 10cm hr-1 given in the Sarmiento&Gruber's book but depending on wind speed 

%================
% DEFINE YOUR SYSTEM
%================

dx=zeros(16,1);
%PO4
dx(1)= -T/V_l*x(1)+T/V_l*x(3)-T/V_l*x(3)*KE_l;
dx(2)= T/V_h*x(1)+M/V_h*x(3)-M/V_h*x(2)-T/V_h*x(2)-T/V_h*x(1)*KE_h-M/V_h*x(3)*KE_h;
dx(3)= T/V_d*x(2)+M/V_d*x(2)-M/V_d*x(3)-T/V_d*x(3)+T/V_d*x(3)*KE_l+T/V_d*x(1)*KE_h+M/V_d*x(3)*KE_h;
%DIC
dx(4)=T/V_l*x(6)-T/V_l*x(4)-T/V_l*KE_l*Rcp*x(3)-T/V_l*KE_l*Rcp*x(3)*Rpump -(co2_flux(Kwind, temp_l, salinity_l,x(4), x(8), x(7),wK)*area_l)*(60*60*24*365)/V_l;
dx(5)=T/V_h*x(4)-T/V_h*x(5)+M/V_h*x(6)-M/V_h*x(5)-T/V_h*KE_h*Rcp*x(1)-M/V_h*KE_h*Rcp*x(3)-(co2_flux(Kwind, temp_h, salinity_h,x(5), x(9), x(7),wK)*area_h)*(60*60*24*365)/V_h;
dx(6)=T/V_d*x(5)-T/V_d*x(6)+M/V_d*x(5)-M/V_d*x(6)+T/V_d*KE_l*Rcp*x(3)+T/V_d*KE_h*Rcp*x(1)+M/V_d*KE_h*Rcp*x(3)+T/V_d*KE_l*Rcp*x(3)*Rpump;
if isnan(setCO2)
	%disp('Keeping ddiagnostic CO2')
	dx(7)=((co2_flux(Kwind, temp_l, salinity_l,x(4), x(8), x(7),wK))*area_l)*(60*60*24*365)/V_a + ((co2_flux(Kwind, temp_h, salinity_h, x(5), x(9), x(7),wK))*area_h)*(60*60*24*365)/V_a;
 else
 	%disp('Holding CO2 levels constant at initialised initial condition')
	dx(7)=0;
end
% alkalinity
dx(8)= T/V_l*x(10)-T/V_l*x(8)+T/V_l*KE_l*Rnp*x(3)-T/V_l*KE_l*Rcp*x(3)*Rpump*Ralk; 
dx(9)= T/V_h*x(8)+M/V_h*x(10)-M/V_h*x(9)-T/V_h*x(9)+T/V_h*KE_h*Rnp*x(1)+M/V_h*KE_h*Rnp*x(3); 
dx(10)= T/V_d*x(9)+M/V_d*x(9)-M/V_d*x(10)-T/V_d*x(10)-T/V_d*KE_l*Rnp*x(3)-T/V_d*KE_h*Rnp*x(1)-M/V_d*KE_h*Rnp*x(3)+T/V_d*KE_l*Rcp*x(3)*Rpump*Ralk; 
%temperature
dx(11)=0; % set LL surf temp
dx(12)=0; % set HL surf temp
dx(13)= T/V_d*x(12)+M/V_d*x(12)-M/V_d*x(13)-T/V_d*x(13);
%salinity
dx(14)=0; % set LL surf Sal
dx(15)=0; % set LL surf Sal
dx(16)=T/V_d*x(15)+M/V_d*x(15)-M/V_d*x(16)-T/V_d*x(16);
% CaCO3 compensation
if isnan(setSScsh)
	%disp('No CaCO3 compensation; closed-system')	
 else
 	%disp('Holding CO2 levels constant at initialised initial condition')
	CSH = carb_solver(x(13)+273.15,x(16),x(6), x(10), setSScsh,wK) %deep ocean carbonate chemistry
	netCaCO3dissolution = -(CSH-setSScsh) * (0.02/100) * (40e12)/V_d
	dx(6) = dx(6)   + 1*netCaCO3dissolution/V_d; %DIC from net CaCO3 dissolution/preservation
	dx(10) = dx(10) + 2*netCaCO3dissolution/V_d; %ALK from net CaCO3 dissolution/preservation
end
end 



