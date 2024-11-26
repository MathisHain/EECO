KE_h = 0.1;
ALKmean = 2364*10^-6;%(mol/kg)
DICmean = 2249*10^-6;%(mol/kg)
whichK='modern';
setSScsh = 2657; %set the desired steady state CSH depth for CaCO3 compensation; if NaN its closed-system CaCO3
setCO2 = NaN; % spin-up to a certain CO2; if NaN its closed-system CO2
finalstate = boxmodel4_function(KE_h,ALKmean,DICmean,whichK,setSScsh,setCO2);

N = length(finalstate(:,8));
CSH = zeros(N,1);
for id = 1:N
    CSH(id) = carb_solver(finalstate(id,13),finalstate(id,16),finalstate(id,6), finalstate(id,10), 3000,whichK);
end

% x vector:
% 1 = PO4_ll; 2 = PO4_hl; 3 = PO4_d; 4 = DIC_ll; 5 =  DIC_hl; 6 = DIC_D; 
% 7 = pCO2_a; 8 = Alk_ll; 9 = ALk_hl; 10 = Alk_d; 11 = T_ll; 12 = T_hl
% 13 = T_d; 14 = S_ll; 15 = S_hl; 16 = S_d

fprintf('Final CO2: %d\n',finalstate(end,7))
fprintf('Final H-PO4: %d\n',finalstate(end,2))
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
plot(CSH) % CSH

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
