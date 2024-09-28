KE_h = 0.2;
ALKmean = 2364*10^-6;%(mol/kg)
DICmean = 2249*10^-6;%(mol/kg)
whichK='modern';
setSScsh = NaN; %set the desired steady state CSH depth for CaCO3 compensation; if NaN its closed-system CaCO3
setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2);

N = length(finalstate(:,8));
CSH = zeros(N,1);
for id = 1:N
    CSH(id) = carb_solver(finalstate(id,13)+273.15,finalstate(id,16),finalstate(id,6), finalstate(id,10), 4000,whichK);
end

% x vector:
% 1 = PO4_ll; 2 = PO4_hl; 3 = PO4_d; 4 = DIC_ll; 5 =  DIC_hl; 6 = DIC_D; 
% 7 = pCO2_a; 8 = Alk_ll; 9 = ALk_hl; 10 = Alk_d; 11 = T_ll; 12 = T_hl
% 13 = T_d; 14 = S_ll; 15 = S_hl; 16 = S_d

fprintf('Final CO2: %d\n',finalstate(end,7))
fprintf('Final deepDIC: %d\n',finalstate(end,6))
fprintf('Final deepALK: %d\n',finalstate(end,10))
fprintf('Final CSH: %d\n',carb_solver(finalstate(end,13)+273.15,finalstate(end,16),finalstate(end,6), finalstate(end,10), 4000,whichK))


% now a loop to run for different levels of Npaz for different combinations of whichK and setCSH

h=figure;
subplot(4,1,1)
plot(finalstate(:,7)) %CO2

subplot(4,1,2)
plot(CSH)

subplot(4,1,3)
hold on
plot(finalstate(:,6)) % deep DIC
plot(finalstate(:,10)) % deep ALK
hold off 

subplot(4,1,4)
hold on
plot(finalstate(:,11)) % LL T
plot(finalstate(:,12)) % HL T
plot(finalstate(:,13)) % deep T
hold off 

saveas(h,'QuickPlot','jpg')
