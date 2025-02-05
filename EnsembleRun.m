clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START MODERN ENSEMBLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Running MODERN open-system spin-up\n')
% Determine deltaALKmean for modern open-system reference case
KE_h = 0.1;
ALKmean = 2364*10^-6;%(mol/kg)
DICmean = 2255*10^-6;%(mol/kg)
whichK='modern';
setSScsh = 3000; %set the desired steady state CSH depth for CaCO3 compensation; if NaN its closed-system CaCO3
setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
Tfeedback = 1;
init_dT = 0;
spinupM = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,50000);
deltaALK=spinupM(end,18)-ALKmean;
init_dT = spinupM(end,11)-273.15-25;
fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',spinupM(end,18)-ALKmean , spinupM(end,7) , spinupM(end,17),spinupM(end,11)-273.15)

modernK_closed = {};
for polarPid = 1:1
	KE_h = 0.1*polarPid;
	fprintf('Running MODERN/CLOSED experiment with KE_h=%d\n',KE_h)
	ALKmean = 2364*10^-6 + deltaALK;%(mol/kg)
	DICmean = 2255*10^-6 + deltaALK/2;%(mol/kg)
	setSScsh = NaN; %set the desired steady state CSH depth for CaCO3 compensation; if NaN its closed-system CaCO3
	setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
	Tfeedback = 0;
	finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,10000);
	fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',finalstate(end,18)-ALKmean , finalstate(end,7) , finalstate(end,17),finalstate(end,11)-273.15)
	modernK_closed{end+1} = finalstate;
end

modernK_open = {};
for polarPid = 1:1
	KE_h = 0.1*polarPid;
	fprintf('Running MODERN/OPEN experiment with KE_h=%d\n',KE_h)
	ALKmean = 2364*10^-6 + deltaALK;%(mol/kg)
	DICmean = 2255*10^-6 + deltaALK/2;%(mol/kg)
	setSScsh = 3000;
	setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
	Tfeedback = 0;
	finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,10000);
	fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',finalstate(end,18)-ALKmean , finalstate(end,7) , finalstate(end,17),finalstate(end,11)-273.15)
	modernK_open{end+1} = finalstate;
end

modernK_open_Tfeed = {};
for polarPid = 1:1
	KE_h = 0.1*polarPid;
	fprintf('Running MODERN/OPEN/TFEED experiment with KE_h=%d\n',KE_h)
	ALKmean = 2364*10^-6 + deltaALK;%(mol/kg)
	DICmean = 2255*10^-6 + deltaALK/2;%(mol/kg)
	setSScsh = 3000;
	setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
	Tfeedback = 1;
	finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,10000);
	fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',finalstate(end,18)-ALKmean , finalstate(end,7) , finalstate(end,17),finalstate(end,11)-273.15)
	modernK_open_Tfeed{end+1} = finalstate;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MODERN ENSEMBLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START EOCENE ENSEMBLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('Running EOCENE open-system spin-up\n')
% Determine deltaALKmean for modern open-system reference case
KE_h = 0.1;
ALKmean = 2364*10^-6;%(mol/kg)
DICmean = 2255*10^-6;%(mol/kg)
whichK='Eocene';
setSScsh = 3000; %set the desired steady state CSH depth for CaCO3 compensation; if NaN its closed-system CaCO3
setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
Tfeedback = 1;
init_dT = 0;
spinupE = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,50000);
deltaALK=spinupE(end,18)-ALKmean;
init_dT = spinupE(end,11)-273.15-25;
fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',spinupE(end,18)-ALKmean , spinupE(end,7) , spinupE(end,17),spinupE(end,11)-273.15)

eoceneK_closed = {};
for polarPid = 1:1
	KE_h = 0.1*polarPid;
	fprintf('Running EOCENE/CLOSED experiment with KE_h=%d\n',KE_h)
	ALKmean = 2364*10^-6 + deltaALK;%(mol/kg)
	DICmean = 2255*10^-6 + deltaALK/2;%(mol/kg)
	setSScsh = NaN; %set the desired steady state CSH depth for CaCO3 compensation; if NaN its closed-system CaCO3
	setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
	Tfeedback = 0;
	finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,10000);
	fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',finalstate(end,18)-ALKmean , finalstate(end,7) , finalstate(end,17),finalstate(end,11)-273.15)
	modernK_closed{end+1} = finalstate;
end

eoceneK_open = {};
for polarPid = 1:1
	KE_h = 0.1*polarPid;
	fprintf('Running EOCENE/OPEN experiment with KE_h=%d\n',KE_h)
	ALKmean = 2364*10^-6 + deltaALK;%(mol/kg)
	DICmean = 2255*10^-6 + deltaALK/2;%(mol/kg)
	setSScsh = 3000;
	setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
	Tfeedback = 0;
	finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,10000);
	fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',finalstate(end,18)-ALKmean , finalstate(end,7) , finalstate(end,17),finalstate(end,11)-273.15)
	eoceneK_open{end+1} = finalstate;
end

eoceneK_open_Tfeed = {};
for polarPid = 0:2
	KE_h = 0.1*polarPid;
	fprintf('Running EOCENE/OPEN/TFEED experiment with KE_h=%d\n',KE_h)
	ALKmean = 2364*10^-6 + deltaALK;%(mol/kg)
	DICmean = 2255*10^-6 + deltaALK/2;%(mol/kg)
	setSScsh = 3000;
	setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
	Tfeedback = 1;
	finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,10000);
	fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',finalstate(end,18)-ALKmean , finalstate(end,7) , finalstate(end,17),finalstate(end,11)-273.15)
	eoceneK_open_Tfeed{end+1} = finalstate;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START EOCENE ENSEMBLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x vector:
% 1 = PO4_ll; 2 = PO4_hl; 3 = PO4_d; 4 = DIC_ll; 5 =  DIC_hl; 6 = DIC_D; 
% 7 = pCO2_a; 8 = Alk_ll; 9 = ALk_hl; 10 = Alk_d; 11 = T_ll; 12 = T_hl
% 13 = T_d; 14 = S_ll; 15 = S_hl; 16 = S_d
finalstate = spinupE;

fprintf('Final CO2: %d\n',finalstate(end,7))
fprintf('Final L-PO4: %d\n',finalstate(end,1))
fprintf('Final H-PO4: %d\n',finalstate(end,2))
fprintf('Final D-PO4: %d\n',finalstate(end,3))
fprintf('Final deepDIC: %d\n',finalstate(end,6))
fprintf('Final deepALK: %d\n',finalstate(end,10))
fprintf('Final Tk: %d\n',finalstate(end,13))
fprintf('Final Sal: %d\n',finalstate(end,16))
fprintf('Final CSH: %d\n',carb_solver(finalstate(end,13),finalstate(end,16),finalstate(end,6), finalstate(end,10), 3000,whichK))


% now a loop to run for different levels of Npaz for different combinations of whichK and setCSH

h=figure;
subplot(3,2,1)
plot(finalstate(:,7)*1000000) %pCO2

subplot(3,2,3)
hold on
plot(finalstate(:,4)*1000) % LL DIC
plot(finalstate(:,5)*1000) % HL DIC
plot(finalstate(:,6)*1000) % deep DIC
hold off 

subplot(3,2,5)
hold on
plot(finalstate(:,8)*1000) % LL ALK
plot(finalstate(:,9)*1000) % HL ALK
plot(finalstate(:,10)*1000) % deep ALK
hold off 

subplot(3,2,2)
hold on
plot(finalstate(:,1)*1000000) % LL PO4
plot(finalstate(:,2)*1000000) % HL PO4
plot(finalstate(:,3)*1000000) % deep PO4
plot(finalstate(:,17)) % CSH

hold off 

subplot(3,2,4)
hold on
plot(finalstate(:,11)-273.15) % LL T
plot(finalstate(:,12)-273.15) % HL T
plot(finalstate(:,13)-273.15) % deep T
hold off 

subplot(3,2,6)
hold on
plot(finalstate(:,14)) % LL S
plot(finalstate(:,15)) % HL S
plot(finalstate(:,16)) % deep S
hold off 

saveas(h,'QuickPlot','jpg')
