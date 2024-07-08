Npaz = 30*10^-6;%(mol/kg)
ALKmean = 2364*10^-6;%(mol/kg)
DICmean = 2249*10^-6;%(mol/kg)
whichK='modern';
setSScsh = NaN; %set the desired steady state CSH depth for CaCO3 compensation; if NaN its closed-system CaCO3
setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
finalstate = boxmodel4_function(Npaz,ALKmean,DICmean,whichK,setSScsh,setCO2);

N = length(finalstate(:,8));
CSH = zeros(N,1);
for id = 1:N
    CSH(id) = carb_solver(finalstate(id,22)+273.15,finalstate(id,25),finalstate(id,7), finalstate(id,19), 4000,whichK);
end

fprintf('Final CO2: %d\n',finalstate(end,8))
fprintf('Final deepDIC: %d\n',finalstate(end,7))
fprintf('Final deepALK: %d\n',finalstate(end,19))
fprintf('Final CSH: %d\n',carb_solver(finalstate(end,22)+273.15,finalstate(end,25),finalstate(end,7), finalstate(end,19), 4000,whichK))


% now a loop to run for different levels of Npaz for different combinations of whichK and setCSH

h=figure;
subplot(4,1,1)
plot(finalstate(:,8)) %CO2

subplot(4,1,2)
plot(CSH)

subplot(4,1,3)
hold on
plot(finalstate(:,7)) % deep DIC
plot(finalstate(:,19)) % deep ALK
hold off 

subplot(4,1,4)
hold on
plot(finalstate(:,20)) % LL T
plot(finalstate(:,21)) % HL T
plot(finalstate(:,22)) % deep T
hold off 

saveas(h,'QuickPlot','jpg')