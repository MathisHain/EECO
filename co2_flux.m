function [ phico2 ] = co2_flux(Kwind, T,S,DIC, ALK, pco2_a,wK)

%Calculate the pco2 flux between the atmosphere and the surface ocean from
%Temperature, salinity, DIC, Alkalinity of the ocean and the pco2 of the
%air.

%phico2[mole/m2*s] - Carbon flux between the ocean and the atmosphere (output)

%Kwind[m/s] - Wind piston velocity
%Temp[K] - Ocean Temperature
%salinity[g/kg] - Ocean Salinity
%DIC[mole/kg] - Ocean Dissolved inorganic carbon
%ALK[mole/kg] - Ocean Alkalinity
%pco2_a[ppm] - Atmopheric pCO2

period = 1; % 1 for today and 0 for Eocene

 if wK == 'Eocene'

    K0 = exp(-60.2409 + 93.4517*(100.0/T)+23.3585*log(T/100.0) + S*(0.023517 - 0.023656*(T/100.0) + 0.0047036*(T/100.0)^2)); %% mol kg-1 atm-1 (Weiss, 1974; table 8.2.2 in Sarmiento&Gruber) Hain 2015

    K1 = 10.0^(61.2172 - 3633.86/T - 9.6777*log(T) + 0.011555*S - 0.0001152*S^2); %% mol kg-1 Hain 2015

    K2 = 10.0^(-25.929 - 471.78/T  + 3.16967*log(T) + 0.01781*S - 0.0001122*S^2); %% mol kg-1 Hain 2015

    Kw = exp(148.9652-13847.26/T-23.6521*log(T) + S^(1/2)*(-5.977+118.67/T+1.0495*log(T))-0.01615*S); %mol kg-1 Hain 2015

    Kb = exp(148.0248+137.1942*S^(1/2)+1.62142*S+1/T*(-8966.9-2890.53*S^(1/2)-77.942*S+1.728*S^1.5-0.0996*S^2)+log(T)*(-24.4344-25.085*S^(1/2)-0.2474*S)+0.053105*S^(1/2)*T); %(mol kg-1)2 %mol kg-1 Hain 2015
    
 elseif wK == 'modern'
    
    K0 = exp(-59.10558297 + 91.85983392*(100.0/T)+22.80815968*log(T/100.0) + S*(0.026034049 - 0.025313942*(T/100.0) + 0.004987136*(T/100.0)^2)); %% mol kg-1 atm-1 (Weiss, 1974; table 8.2.2 in Sarmiento&Gruber) Hain 2015

    K1 = 10.0^(62.19335659 - 3669.526066/T - 9.834540861*log(T) + 0.011193361*S - 0.000113713*S^2); %% mol kg-1 Hain 2015

    K2 = 10.0^(-26.71121142 - 437.8112533/T  + 3.287439691*log(T) + 0.017419883*S - 0.000110534*S^2); %% mol kg-1 Hain 2015

    Kw = exp(146.6165292-13701.93677/T-23.34510121*log(T) + S^(1/2)*(-5.777126033+105.5781367/T+1.023014734*log(T))-0.016256003*S); %mol kg-1 Hain 2015

    Kb = exp(156.169978+138.6678081*S^(1/2)+1.744675879*S+1/T*(-9226.072384-2874.822309*S^(1/2)-82.05134674*S+1.444274205*S^1.5-0.09130297*S^2)+log(T)*(-25.71609713-25.43581085*S^(1/2)-0.264978916*S)+0.054449827*S^(1/2)*T); %(mol kg-1)2 %mol kg-1 Hain 2015
    
 end
 
c = 11.88*10^(-6);%11.88*10^(-6)% not the same as in Sarmiento&Gruber, c = 1.185*10^(-6) (Uppstrom, 1974)%give the relationship between borate and salinity

nt=10000;


PH_firstguess = 6.5;
H = 10^(-PH_firstguess); % mol kg-1 

for i = 1:nt % see table 8.2.1 in Sarmiento&Gruber. To solve H+ using an iterative approach. 

    guessed = DIC/(H/K1+1+K2/H)+2*(DIC/(H^2/(K1*K2)+H/K2+1))+Kw/H-H+Kb*c*S/(H+Kb)  ;  % mol kg-1  to see table 
    
    if (abs(guessed-ALK) < 10^-6)
      break
    elseif (guessed-ALK < 0)
        PH_firstguess = PH_firstguess + 0.001;
        H = 10^(-PH_firstguess);
    elseif (guessed-ALK > 0)
        PH_firstguess = PH_firstguess - 0.001;
        H = 10^(-PH_firstguess);
    end  
    
% Once H+ has been solved, we can calculate pCO2 from DIC and H+ as follow (see
% table 8.2.1 in Sarmiento&Gruber
    
pco2=DIC/K0*(H^2/(H^2+K1*H+K1*K2)); %atm and 280ppm = 280µatm


phico2 = -K0*Kwind*(pco2_a-pco2)*1000; %mol m-2 s-1
end