%%%%%%%     MODEL MAIN TEMPLATE    %%%%%%%
%==============================================
%INITIAL CONDITIONS and time span of simulation
%===============================================
    
PO4_ll_ini = 0.5*10^-6;%(mol/kg)
PO4_hl_ini = 1.48*10^-6;%(mol/kg)
PO4_d_ini  = 2.15*10^-6;%(mol/kg)
PO4_a_ini  = 0.0; 
DIC_ll_ini = 1924*10^-6;%(mol/kg)1924*10^-6%%%%3000*10^-6
DIC_hl_ini = 2149*10^-6;%(mol/kg)2149*10^-6%%%%3351*10^-6
DIC_D_ini  = 2240*10^-6;%(mol/kg)2240*10^-6%%%%3493*10^-6
pCO2_a_ini = 280*10^-6 ;%(atm)) 
NO3_ll_ini = 5*10^-6;%(mol/kg)
NO3_hl_ini = 25*10^-6;%(mol/kg)
NO3_d_ini  = 34.7*10^-6;%(mol/kg)
NO3_a_ini  = 0.0;   
global ALK_l ALK_h
ALK_l   = 2322*10^-6; %Low latitude alkalinity [mol/kg] - 2322*10^-6
ALK_h   = 2322*10^-6; %high latitude alkalinity [mol/kg] - 2322*10^-6

d15N_ll_ini        = 5.4; % d15N for delta values
d15N_hl_ini        = 5.4;
d15N_d_ini         = 5.4;
d15N_a_ini         = 5.4;
R15ref             = 0.0036782; % isotopic ratio of Air N2 (i.e., the international reference for N isotopes)

N15_ll_ini         =((((d15N_ll_ini)/1000)+1)*R15ref)*NO3_ll_ini; %15N for 15N concentration
N15_hl_ini         =((((d15N_hl_ini)/1000)+1)*R15ref)*NO3_hl_ini;
N15_d_ini          =((((d15N_d_ini)/1000)+1)*R15ref)*NO3_d_ini;
N15_a_ini          =((((d15N_a_ini)/1000)+1)*R15ref)*NO3_a_ini;

x0 = [PO4_ll_ini, PO4_hl_ini, PO4_d_ini, PO4_a_ini, DIC_ll_ini, DIC_hl_ini,DIC_D_ini,pCO2_a_ini,NO3_ll_ini,NO3_hl_ini,NO3_d_ini,NO3_a_ini,N15_ll_ini,N15_hl_ini,N15_d_ini,N15_a_ini];

tspan = (0:1:5000); %1000 years of simulation
%=================
%SOLVING THE ODE
%=================


    [t,x] = ode45(@CO2atm_ode,tspan,x0,[]);
    d15N=zeros(length(tspan),3);
    d15N(:,1)     = (((x(:,13)./x(:,9))./R15ref)-1)*1000;
    d15N(:,2)     = (((x(:,14)./x(:,10))./R15ref)-1)*1000;
    d15N(:,3)     = (((x(:,15)./x(:,11))./R15ref)-1)*1000;
    

%==============
%plot results
%==============
t_yr=t; %years in second/n. of second in a year = years
%plot for phosphorous
%figure(1);
%hold on;
%plot (t_yr,x(:,1:3));
%title('Phosphorous vs time');
%xlabel('time (years)');
%ylabel ('mol PO4/Kg');
%legend ('PO4_l','PO4_h','PO4_d');
%hold on;
%plot for nitrogen
%figure(2);
%hold on;
%plot (t_yr,x(:,9:11));
%title('nitrogen vs time');
%xlabel('time (years)');
%ylabel ('mol NO3/Kg');
%legend ('NO3_l','NO3_h','NO3_d');
%hold on;
%plot for d15N
%figure(3);
%hold on;
%plot (t_yr,d15N(:,1:3));
%title('d15N vs time');
%xlabel('time (years)');
%ylabel ('d15N');
%legend ('d15N_l','d15N_h','d15N_d');
%hold on;
%plot for carbon ------maybe change time from seconds to years?
figure(4);
hold on;
plot(t_yr,x(:,5:7));
title('DIC vs time');
xlabel('time (years)');
ylabel ('mol C/Kg');
legend('DIC_l','DIC_h','DIC_D');
hold on;
%plot for pCO2
%x8ppm = x(8)/((1027/1000000)*(8.205736*10^-5)*273)/10e-6
figure(5);
hold on;
plot(t_yr,x(:,8)*10^6);
title ('CO2 vs time');
xlabel ('time (years)');
ylabel ('CO2 ppmv');
hold on;

DIC_LL = x(length(tspan),5).*1000000;
DIC_HL = x(length(tspan),6).*1000000;
PCO2ATM = x(length(tspan),8).*1000000;

%to give alkalinity
[A] = CO2SYS(ALK_l*1000000,DIC_LL,1,2,34,0,15,4200,0,15,1,1,5,1);%
PH_LL = A(1,18);
[B] = CO2SYS(ALK_h*1000000,DIC_HL,1,2,31,0,2,4200,0,15,1,1,5,1);%
PH_HL = B(1,18);

clearvars -except A B DIC_LL DIC_HL PCO2ATM PH_LL PH_HL


