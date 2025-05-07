
clear
%load('devrun')
if exist('spinupM','var')
	fprintf('REUSING OUTPUT DATA ==> SAVES TIME\n')
else
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START MODERN ENSEMBLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ex = 0:10;
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
	fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',spinupM(end,18)-ALKmean, spinupM(end,7), spinupM(end,17),spinupM(end,11)-273.15)

	modernK_closed = {};
	for polarPid = ex
		KE_h = 0.1*polarPid;
		fprintf('Running MODERN/CLOSED experiment with KE_h=%d\n',KE_h)
		ALKmean = 2364*10^-6 + deltaALK;%(mol/kg)
		DICmean = 2255*10^-6 + deltaALK/2;%(mol/kg)
		setSScsh = NaN; %set the desired steady state CSH depth for CaCO3 compensation; if NaN its closed-system CaCO3
		setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
		Tfeedback = 0;
		finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,20000);
		fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',finalstate(end,18)-ALKmean , finalstate(end,7) , finalstate(end,17),finalstate(end,11)-273.15)
		modernK_closed{end+1} = finalstate;
	end

	modernK_closed_Tfeed= {};
	for polarPid = ex
		KE_h = 0.1*polarPid;
		fprintf('Running MODERN/OPEN experiment with KE_h=%d\n',KE_h)
		ALKmean = 2364*10^-6 + deltaALK;%(mol/kg)
		DICmean = 2255*10^-6 + deltaALK/2;%(mol/kg)
		setSScsh = NaN;
		setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
		Tfeedback = 1;
		finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,20000);
		fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',finalstate(end,18)-ALKmean , finalstate(end,7) , finalstate(end,17),finalstate(end,11)-273.15)
		modernK_closed_Tfeed{end+1} = finalstate;
	end

	modernK_open_Tfeed = {};
	for polarPid = ex
		KE_h = 0.1*polarPid;
		fprintf('Running MODERN/OPEN/TFEED experiment with KE_h=%d\n',KE_h)
		ALKmean = 2364*10^-6 + deltaALK;%(mol/kg)
		DICmean = 2255*10^-6 + deltaALK/2;%(mol/kg)
		setSScsh = 3000;
		setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
		Tfeedback = 1;
		finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2,Tfeedback,init_dT,20000);
		fprintf('∆ALK= %d, CO2=%d, CSH=%d, T=%d\n\n',finalstate(end,18)-ALKmean , finalstate(end,7) , finalstate(end,17),finalstate(end,11)-273.15)
		modernK_open_Tfeed{end+1} = finalstate;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MODERN ENSEMBLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('devrun')


% x vector:
% 1 = PO4_ll; 2 = PO4_hl; 3 = PO4_d; 4 = DIC_ll; 5 =  DIC_hl; 6 = DIC_D; 
% 7 = pCO2_a; 8 = Alk_ll; 9 = ALk_hl; 10 = Alk_d; 11 = T_ll; 12 = T_hl
% 13 = T_d; 14 = S_ll; 15 = S_hl; 16 = S_d
% 17 = CSH
% 18 = ALKmean

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COLLECT STEADY STATE OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(modernK_closed)
MCC = zeros(N,6);
MCT = zeros(N,6);
MOT = zeros(N,6);

for n = 1:N
	MCC(n,1) = modernK_closed{n}(end,7).*10^6; 			% CO2
	MCC(n,2) = modernK_closed{n}(end,11)-273.15;	% Tll
	MCC(n,3) = modernK_closed{n}(end,17);			% CSH
	MCC(n,4) = modernK_closed{n}(end,18);			% ALKmean
	MCC(n,5) = modernK_closed{n}(end,2).*10^6;			% PO4hl
	MCC(n,6) = modernK_closed{n}(end,3).*10^6;			% PO4d
	
	MCT(n,1) = modernK_closed_Tfeed{n}(end,7).*10^6; 			% CO2
	MCT(n,2) = modernK_closed_Tfeed{n}(end,11)-273.15;	% Tll
	MCT(n,3) = modernK_closed_Tfeed{n}(end,17);			% CSH
	MCT(n,4) = modernK_closed_Tfeed{n}(end,18);			% ALKmean
	MCT(n,5) = modernK_closed_Tfeed{n}(end,2).*10^6;			% PO4hl
	MCT(n,6) = modernK_closed_Tfeed{n}(end,3).*10^6;			% PO4d	
	
	MOT(n,1) = modernK_open_Tfeed{n}(end,7).*10^6; 			% CO2
	MOT(n,2) = modernK_open_Tfeed{n}(end,11)-273.15;	% Tll
	MOT(n,3) = modernK_open_Tfeed{n}(end,17);			% CSH
	MOT(n,4) = modernK_open_Tfeed{n}(end,18);			% ALKmean
	MOT(n,5) = modernK_open_Tfeed{n}(end,2).*10^6;			% PO4hl
	MOT(n,6) = modernK_open_Tfeed{n}(end,3).*10^6;			% PO4d		
end
	
	



%%% quickplot
if (1)
	
	
	h = figure;
	for p = 1:4
		subplot(2,2,p)
		hold on
		plot(MCC(:,5),MCC(:,p),"o-",'color','b')
		plot(MCT(:,5),MCT(:,p),"o--",'color','b')
		plot(MOT(:,5),MOT(:,p),"o:",'color','b')

	end
	

	

	saveas(h,'EnsemblePlot','jpg')
end

%%% quick plot
if (0)
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

	h          = figure;
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
end
