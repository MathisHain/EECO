Npaz = 25*10^-6;%(mol/kg)
ALKmean = 2322*10^-6;%(mol/kg)
whichK='modern';
setSScsh = 3500; %set the desired steady state CSH depth for CaCO3 compensation
setCO2 = NaN; % spin-up to a certain CO2
finalstate = boxmodel4_function(Npaz,ALKmean,whichK,setSScsh,setCO2);

finalstate(:,8) %CO2