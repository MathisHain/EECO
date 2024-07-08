function [CSH,H, pH, pCO2, H2CO3, HCO3, CO3, Ozd] = carb_solver(T,S,DIC, ALK, Depth,wK)

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

%period = 1; % 1 for today and 0 for Eocene

if wK == 'modern'
 	Ca= 0.0102821; Mg = 0.0528171;

%{
	 ------------------INTERACTIVE_MODE----------------------
	 --- MyAMI specific ion interaction model Version 0.9 ---
	 ------------------INTERACTIVE_MODE----------------------

	 Enter [Ca2+] in mol/kg for nominal Sal=35 (modern 0.0102821): 0.0102821
	 Enter [Mg2+] in mol/kg for nominal Sal=35 (modern 0.0528171): 0.0528171

	 CODE IS RUNNING ... be patient. Sometimes it takes a while to fit the Temp/Sal function to the Model output.

	 lnK0 = -60.24 (+) 93.45(100/T) (+) 23.36*ln(T/100) (+) Sal*[ 0.023517 (+) -0.023656*(T/100) (+) 0.0047036*(T/100)*(T/100)]
	 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1. 1. 1. 1. 1. 1.]
	 Average accuracy of fitting the function to the model output: 1.0

	 log10K1 = 61.2172 (+) -3633.86/T (+) -9.6777*lnT (+) 0.011555*Sal (+) -0.0001152*Sal*Sal
	 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1. 1. 1. 1. 1.]
	 Average accuracy of fitting the function to the model output: 1.0

	 log10K2 = -25.9290 (+) -471.78/T (+) 3.16967*lnT (+) 0.01781*Sal (+) -0.0001122*Sal*Sal
	 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1. 1. 1. 1. 1.]
	 Average accuracy of fitting the function to the model output: 0.9999999999999998

	 lnKb = 148.0248 (+) 137.1942*Sal^0.5 (+) 1.62142*Sal (+) (1/T*[-8966.90 (+) -2890.53*sqrtSal (+) -77.942*Sal (+) 1.728*Sal^1.5 (+) -0.0996*Sal^2] (+) lnT*[-24.4344 (+) -25.085*Sal^0.5 (+) -0.2474*Sal] (+) 0.053105*T*Sal^0.5
	 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
	 Average accuracy of fitting the function to the model output: 1.0000000000000002

	 lnKw = 148.9652 (+) -13847.26/T (+) -23.6521*lnT (+) Sal^0.5*[118.67/T (+) -5.977 (+) 1.0495*lnT] (+) -0.01615*Sal
	 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1. 1. 1. 1. 1. 1. 1.]
	 Average accuracy of fitting the function to the model output: 1.0

	 log10KspC = -171.9065 (+) -0.077993*T (+) 2839.319/T (+) 71.595*log10T (+) Sal^0.5*(-0.777120 (+) 0.0028426*T (+) 178.340/T) (+) -0.07711*Sal (+) 0.0041249*Sal^1.5
	 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1. 1. 1. 1. 1. 1. 1. 1. 1.]
	 Average accuracy of fitting the function to the model output: 1.0

	 log10KspA = -171.9450 (+) -0.077993*T (+) 2903.293/T (+) 71.595*log10T (+) Sal^0.5*(-0.068393 (+) 0.0017276*T (+) 88.135/T) (+) -0.10018*Sal (+) 0.0059415*Sal^1.5
	 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1. 1. 1. 1. 1. 1. 1. 1. 1.]
	 Average accuracy of fitting the function to the model output: 1.0

	 lnKSO4 = 141.328 (+) -4276.1/T (+) -23.093*lnT (+) (I^0.5*[-13856/T (+) 324.57 (+) -47.986lnT] (+) (I*[35474/T (+) -771.54 (+) 114.723*lnT] (+) -2698/T*I^1.5 (+) 1776/T*I^2 + ln(1-0.001005*Sal)
	 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
	 Average accuracy of fitting the function to the model output: 1.0
%}
	 
    K0 = exp(-60.2409 + 93.4517*(100.0/T)+23.3585*log(T/100.0) + S*(0.023517 - 0.023656*(T/100.0) + 0.0047036*(T/100.0)^2)); %% mol kg-1 atm-1 (Weiss, 1974; table 8.2.2 in Sarmiento&Gruber) Hain 2015

    K1 = 10.0^(61.2172 - 3633.86/T - 9.6777*log(T) + 0.011555*S - 0.0001152*S^2); %% mol kg-1 Hain 2015

    K2 = 10.0^(-25.929 - 471.78/T  + 3.16967*log(T) + 0.01781*S - 0.0001122*S^2); %% mol kg-1 Hain 2015

    Kw = exp(148.9652-13847.26/T-23.6521*log(T) + S^(1/2)*(-5.977+118.67/T+1.0495*log(T))-0.01615*S); %mol kg-1 Hain 2015

    Kb = exp(148.0248+137.1942*S^(1/2)+1.62142*S+1/T*(-8966.9-2890.53*S^(1/2)-77.942*S+1.728*S^1.5-0.0996*S^2)+log(T)*(-24.4344-25.085*S^(1/2)-0.2474*S)+0.053105*S^(1/2)*T); %(mol kg-1)2 %mol kg-1 Hain 2015

	KspC = 10.0^( -171.9065 - 0.077993*T + 2839.319/T + 71.595*log10(T) + S^0.5*(-0.777120 + 0.0028426*T + 178.340/T) - 0.07711*S + 0.0041249*S^1.5);

elseif wK == 'Eocene'
 	Ca= 0.02; Mg = 0.03;
	 
%{
		 ------------------INTERACTIVE_MODE----------------------
		 --- MyAMI specific ion interaction model Version 0.9 ---
		 ------------------INTERACTIVE_MODE----------------------

		 Enter [Ca2+] in mol/kg for nominal Sal=35 (modern 0.0102821): 0.02
		 Enter [Mg2+] in mol/kg for nominal Sal=35 (modern 0.0528171): 0.03

		 CODE IS RUNNING ... be patient. Sometimes it takes a while to fit the Temp/Sal function to the Model output.

		 lnK0 = -59.11 (+) 91.86(100/T) (+) 22.81*ln(T/100) (+) Sal*[ 0.0260343 (+) -0.0253141*(T/100) (+) 0.004987172*(T/100)*(T/100)]
		 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [0.98115298 0.98296523 0.97643848 1.10704364 1.07009414 1.06028819]
		 Average accuracy of fitting the function to the model output: 0.9999993353966867

		 log10K1 = 60.3593 (+) -3596.46/T (+) -9.5483*lnT (+) 0.011435*Sal (+) -0.0001151*Sal*Sal
		 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [0.98598629 0.98970818 0.98663034 0.9896146  0.99892271]
		 Average accuracy of fitting the function to the model output: 0.9999978492507375

		 log10K2 = -25.3162 (+) -492.86/T (+) 3.06644*lnT (+) 0.01777*Sal (+) -0.0001119*Sal*Sal
		 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [0.97636732 1.0446859  0.96743201 0.99789757 0.99764443]
		 Average accuracy of fitting the function to the model output: 1.0000069324089367

		 lnKb = 151.0931 (+) 134.1726*Sal^0.5 (+) 1.51508*Sal (+) (1/T*[-9122.27 (+) -2829.91*sqrtSal (+) -74.679*Sal (+) 1.845*Sal^1.5 (+) -0.1045*Sal^2] (+) lnT*[-24.8837 (+) -24.520*Sal^0.5 (+) -0.2313*Sal] (+) 0.051792*T*Sal^0.5
		 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1.02072857 0.9779755  0.93441706 1.01732683 0.97902845 0.95813133
		  1.06766533 1.0490566  1.0183897  0.97747899 0.9348774  0.97527975]
		 Average accuracy of fitting the function to the model output: 0.9999999965630919

		 lnKw = 151.9032 (+) -13838.18/T (+) -24.1930*lnT (+) Sal^0.5*[177.04/T (+) -7.478 (+) 1.2750*lnT] (+) -0.01765*Sal
		 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1.01972259 0.99934397 1.02286747 1.49187677 1.25117884 1.21490076
		  1.09268578]
		 Average accuracy of fitting the function to the model output: 1.0000527442529463

		 log10KspC = -64.4350 (+) -0.044654*T (+) 149.620/T (+) 27.778*log10T (+) Sal^0.5*(-0.650175 (+) 0.0026312*T (+) 158.390/T) (+) -0.07707*Sal (+) 0.0041206*Sal^1.5
		 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [0.37482589 0.57254127 0.05269567 0.38799471 0.83664723 0.92563886
		  0.88813442 0.99953527 0.99896835]
		 Average accuracy of fitting the function to the model output: 1.0000000600743981

		 log10KspA = -62.6191 (+) -0.044059*T (+) 168.013/T (+) 27.019*log10T (+) Sal^0.5*(0.061723 (+) 0.0015113*T (+) 67.767/T) (+) -0.10020*Sal (+) 0.0059402*Sal^1.5
		 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [ 0.36418072  0.56491531  0.05786976  0.37738552 -0.90248159  0.874789
		   0.76890325  1.0001543   0.99977839]
		 Average accuracy of fitting the function to the model output: 1.0000001278250659

		 lnKSO4 = 143.391 (+) -4357.0/T (+) -23.404*lnT (+) (I^0.5*[-13942/T (+) 326.07 (+) -48.211lnT] (+) (I*[35387/T (+) -771.08 (+) 114.676*lnT] (+) -2649/T*I^1.5 (+) 1764/T*I^2 + ln(1-0.001005*Sal)
		 Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater: [1.01459377 1.01891997 1.01345344 1.00621594 1.00462797 1.00468338
		  0.99753549 0.99940434 0.99958745 0.98187835 0.99321628]
		 Average accuracy of fitting the function to the model output: 1.000000111918892
%}	
    K0 = exp(-59.10558297 + 91.85983392*(100.0/T)+22.80815968*log(T/100.0) + S*(0.026034049 - 0.025313942*(T/100.0) + 0.004987136*(T/100.0)^2)); %% mol kg-1 atm-1 (Weiss, 1974; table 8.2.2 in Sarmiento&Gruber) Hain 2015

    K1 = 10.0^(62.19335659 - 3669.526066/T - 9.834540861*log(T) + 0.011193361*S - 0.000113713*S^2); %% mol kg-1 Hain 2015

    K2 = 10.0^(-26.71121142 - 437.8112533/T  + 3.287439691*log(T) + 0.017419883*S - 0.000110534*S^2); %% mol kg-1 Hain 2015

    Kw = exp(146.6165292-13701.93677/T-23.34510121*log(T) + S^(1/2)*(-5.777126033+105.5781367/T+1.023014734*log(T))-0.016256003*S); %mol kg-1 Hain 2015

    Kb = exp(156.169978+138.6678081*S^(1/2) + 1.744675879*S+1/T*(-9226.072384-2874.822309*S^(1/2) - 82.05134674*S+1.444274205*S^1.5 - 0.09130297*S^2)+log(T)*(-25.71609713-25.43581085*S^(1/2)-0.264978916*S)+0.054449827*S^(1/2)*T); %(mol kg-1)2 %mol kg-1 Hain 2015
	KspC = 10.0^(-64.4350 - 0.044654*T + 149.620/T + 27.778*log10(T) + S^0.5*(-0.650175 + 0.0026312*T + 158.390/T) - 0.07707*S + 0.0041206*S^1.5);
	
 end
 
 P = 0.1007*Depth; R = 83.144621; % Gas constant; corrected relative to Millero (1995GCA)
 Tke = T; Tc = T-273.15;
 K1d = K1 * exp(-(-25.50 + 0.1271*Tc)/R/Tke*P + 0.5*(-0.00308 + 0.877e-4*Tc)/R/Tke*P*P);
 K2d = K2 * exp(-(-15.82 - 0.0219*Tc)/R/Tke*P + 0.5*(0.00113 - 1.475e-4*Tc)/R/Tke*P*P);
 Kbd = Kb * exp(-(-29.48 + 0.1622*Tc - 0.002608*Tc*Tc)/R/Tke*P + 0.5*(-0.00284)/R/Tke*P*P);
 Kwd = Kw * exp(-(-20.02 + 0.1119*Tc - 0.001409*Tc*Tc)/R/Tke*P + 0.5*(-0.00513 + 0.794e-4*Tc)/R/Tke*P*P);
 %KsAd = KsA * exp(-(-46 + 0.5304*Tc)/R/Tke*P + 0.5*(-0.01176 - 3.692e-4*Tc)/R/Tke*P*P);
 KspCd = KspC * exp(-(-48.76 + 0.5304*Tc)/R/Tke*P + 0.5*(-0.01176 - 3.692e-4*Tc)/R/Tke*P*P);
 
c = 11.88*10^(-6);%11.88*10^(-6)% not the same as in Sarmiento&Gruber, c = 1.185*10^(-6) (Uppstrom, 1974)%give the relationship between borate and salinity

nt=10000;



PH_firstguess = 8;
H = 10^(-PH_firstguess); % mol kg-1 

for i = 1:nt % see table 8.2.1 in Sarmiento&Gruber. To solve H+ using an iterative approach. 

    guessed = DIC/(H/K1d+1+K2d/H)+2*(DIC/(H^2/(K1d*K2d)+H/K2d+1))+Kwd/H-H+Kbd*c*S/(H+Kbd)  ;  % mol kg-1  to see table 
    
    if (abs(guessed-ALK) < 10^-7)
		pH = PH_firstguess;
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
% dO/O = dKsCd/KsCd/dz * dCSH = (+1.8%/100m) *dCSH
% DlogO = DlogKsCd/dz * DCSH = (+1.8%/100m) * DCSH
dlogKsp_dz = 0.018/100; % (+1.8%/100m)
H2CO3=DIC*(H^2/(H^2+K1*H+K1d*K2d));
pCO2=H2CO3/K0; %atm and 280ppm = 280µatm
HCO3 = DIC*(K1d*H/(H^2+K1d*H+K1d*K2d));
CO3 = DIC*(K1d*K2d/(H^2+K1d*H+K1d*K2d));
Ozd = Ca*CO3/KspCd;
CSH = Depth + log(Ozd)/dlogKsp_dz;

end