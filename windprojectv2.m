%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% WIND FARM %%%%%%%
%%% NEW AND RENEWABLE %%%
%%%%%%% ENERGIES %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CODE BY CARLOS PEREZ %%%

%%% Function Windrose 
%%% by Daniel Pereira

clear; close all; clc;

%  ________________________________________________________________________
%% DATA AND PARAMETERS 
% *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *

N_wt = 8; % Number of Wind turbines

PR_wt = 1.6; % Rated power of wind turbines (MW)

Cost = 600; % Wind farm cost (eur/kW)

H = 8.1; %Hub height

rho_air = 1.225; % Air density kg/m3

V_tip = 70; % Tip speed (m/s)

% performance = xlsread('performance.xls','A3:DL25');

% Files with wind data
filename_MERRA = "MERRA.dat";
filename_MAST  = "FinalData.csv";

% Load the MERRA data
wind_MERRA = readtable(filename_MERRA, "Delimiter", "tab");
wind_MERRA = table2array(wind_MERRA);
wind_MERRA = [wind_MERRA(:,6), wind_MERRA(:,5)];

% Load the mast data
wind_MAST = load_gwynt_y_mor(filename_MAST);
wind_MAST = [wind_MAST(:,5), wind_MAST(:,7), ...
    wind_MAST(:,1), wind_MAST(:,2), wind_MAST(:,3), wind_MAST(:,4)];
wind_MAST2 = wind_MAST; NM = size(wind_MAST2, 1); NC = size(wind_MAST2, 2);
wind_MAST_aux = zeros(NM,NC);
% Correct the data from the meteorological mast 
counter = 1;
initial_date = wind_MAST(1, 3:end);
last_date    = wind_MAST(end, 3:end);

year_diff    = last_date(1) - initial_date(1);
day_diff     = last_date(3) - initial_date(3);
month_diff   = last_date(2) - initial_date(2);
if month_diff < 0; year_diff = year_diff - 1; end

limit_date = initial_date; limit_date(1) = initial_date(1) + year_diff;

for ii = 1:NM
    if wind_MAST2(ii, 3:end) == limit_date
        break
    elseif wind_MAST2(ii,1) >= 0
        wind_MAST_aux(counter,1) = wind_MAST2(ii,1);
        wind_MAST_aux(counter,2) = wind_MAST2(ii,2);
        wind_MAST_aux(counter,3) = wind_MAST2(ii,3);
        wind_MAST_aux(counter,4) = wind_MAST2(ii,4);
        wind_MAST_aux(counter,5) = wind_MAST2(ii,5);
        wind_MAST_aux(counter,6) = wind_MAST2(ii,6);
        counter = counter + 1;
    end
end
clear wind_MAST;
wind_MAST = wind_MAST_aux(1:counter-1,1:end);

% Process the data for MERRA
figure_names.weibull = 'weibull_MERRA'; 
figure_names.windrose = 'windrose_MERRA';
wind_analysis(wind_MERRA, figure_names, 23);

% Process the data for the metereological mast
figure_names.weibull = 'weibull_MAST'; 
figure_names.windrose = 'windrose_MAST';
wind_analysis(wind_MAST, figure_names, 37);

return

% ANNUAL AVERAGE WIND POWER DENSITY

v_0 = 0; % Starting speed of the iteration
v_f = 50; % Last speed of the iteration

N = 10000;
delta_v = (v_f - v_0) / N;

PD = zeros(N,1);

speed = 0;
PD(1) = 0;
sumPD = 0;
for i=1:N-1
    speed = speed + delta_v;
    PD(i+1) = k / c .* (speed / c) .^(k-1) .* exp(-(speed/c).^k) * ...
        rho_air * speed^3 / 2;
    sumPD = sumPD + PD(i+1);
end

WPD_numerical = delta_v * ( (PD(1)+PD(N))/2 +sumPD ); 

%% Alternative method for wind power density

WP = zeros(sz_wind(1),1);
WP_mean = 0;
for i=1:sz_wind(1)
    WP(i) = 1/2 * rho_air * v(i,1)^3;
    WP_mean = WP_mean + WP(i);
end

WPD_meanvalues = WP_mean / sz_wind(1);

%% Cp - Lambda curves

v_ci = 4; %cut in velocity m/s
v_co = 25; % cut out velocity m/s

cp_max = 0;

figure;

for j=2:8 % beta iteration
    
    n_values = 23;

    %%% GETTING THE CP MAX FOR EACH VALUE OF BETA
    for ii=1:23 % lambda iteration
        
        if isnan(performance(ii,3+4*(j-1)))
            n_values = n_values - 1;
        elseif performance(ii,3+4*(j-1))>cp_max
            beta_cpmax = performance(1,1+4*(j-1));
            column_lambdaopt = 2+4*(j-1);
            column_cpmax = 3+4*(j-1);
            column_betamax = 1 + 4*(j-1);
            cp_max = performance(ii,3+4*(j-1)); 
            lambda_cpmax = performance(ii,2+4*(j-1));
        end        
    end
   
    plot(performance(1:n_values,2+4*(j-1)),performance(1:n_values,3 +4*(j-1)),'DisplayName',string(performance(2,1+4*(j-1))));
    hold on;

end 
legend;
grid on;
title('Cp - \lambda curves');
xlabel('\lambda');
ylabel('Cp');
        
%% DIAMETER CALCULATION        
  
R_i = 30; %starting radius for iteration
R_f = 35; % last radius for iteration

delta_R = 0.5; % radius step for iteration

N_r = (R_f - R_i)/delta_R + 1; % number of radius for iteration 
 
    N = 10000; % DISCRETIZATION OF VELOCITIES
    delta_v = (v_co - v_ci) / N;

r = R_i; 
   
P_wt = zeros(1,N_r);
D_i = zeros(1,N_r);
CF = zeros(2,N_r);
  
 for iter = 1:N_r 
  
    %%% Calculo de la potencia para Cp max
    speed = v_ci;
    sumPI = 0;
    PI = zeros(N+1,1);
    
    V_R = (PR_wt * 10^6 * 2 / (cp_max * rho_air * pi * r^2))^(1/3); 
    Power = zeros(1,N);
    for i=2:N
        
        if speed < V_R   
            Power(i) = PowerInstantaneous(speed, cp_max, r, rho_air);
        else
            Power(i) = PR_wt*10^6;
        end
        
        f_w = weibull(speed);
        PI(i) = f_w * Power(i);
        
        
        sumPI = sumPI + PI(i);
        speed = speed + delta_v;
    end
    
    Power(1) = PowerInstantaneous(v_ci, cp_max, r, rho_air);
    Power(N+1) = PowerInstantaneous(v_co, cp_max, r, rho_air);
    PI(1) = Power(1) * weibull(v_ci);
    PI(N+1) = Power(N+1) * weibull(v_co);
    P_wt(iter) = delta_v * ( (PI(1)+PI(N+1))/2 + sumPI );
    CF(1, iter) = P_wt(iter)/(PR_wt*10^6);
    CF(2, iter) = r;
    D_i(iter) = 2 * r;
    
    r = r + delta_R;
    
 end
figure;
plot(CF(2,:),CF(1,:),'-o');
title('CF for different values of R');
grid on;
xlabel('R');
ylabel('CF');
 
 %% PLOT OF WIND POWER CURVES
 
 R_wt = 32;
 A_rotor = pi * R_wt^2;
 figure;

 delta_v =  1;
 
 v_range2 = v_ci:delta_v:v_co;
 
 % Wind power
 
 windpower = 1/2 * rho_air * pi * R_wt ^2 * v_range2 .^ 3;
 plot(v_range2,windpower,'DisplayName','Wind Power');
 hold on;
 
 % Betz limit
 
 windpower_betz = 16 / 27 * 1/2 * rho_air * pi * R_wt ^2 * v_range2 .^ 3;
 plot(v_range2,windpower_betz,'DisplayName','Betz Limit');
 xlabel('Wind speed (m/s)');
 ylabel('Power (W)');
 title('Wind Power Curves');
 grid on;
 legend;
 
 %% Constant speed - Stall regulated
 
 delta_v = 1;
 N = (v_co-v_ci)/delta_v + 1;
 
 speed = v_ci;
 
 lambda_stall = zeros(1,N);
 cp_stall = zeros(1,N);
 Power_stall = zeros(1,N);
 Power_stall_rated = zeros(1,N);
 Power_stall_rated_2 = zeros(1,N);
 
 for i = 1:N

     lambda_stall(i) = V_tip / speed;
     
     cp_stall(i) = interp1(performance(1:22,14),performance(1:22,15),lambda_stall(i),'spline');
     
     Power_stall(i) = 1 / 2 * cp_stall(i) * rho_air * pi * R_wt^2 * speed ^3;
     Power_stall_rated(i) = Power_stall(i);
     Power_stall_rated_2(i) = Power_stall(i);
     if Power_stall_rated(i) > PR_wt * 10^6
         Power_stall_rated(i) = PR_wt * 10^6;
         Power_stall_rated_2(i) = 0;
     end
     
    speed = speed + delta_v; 
 end
 
 figure; 
 plot(v_ci:delta_v:v_co,Power_stall_rated,'DisplayName','Rated power');
 title('Power output of stall wind turbine')
 grid on;
 xlabel('Wind speed (m/s)');
 ylabel('Power (W)');
 hold on;
 plot(v_ci:delta_v:v_co,Power_stall,'DisplayName','Stall wind turbine');
 legend;
 
 figure;
 subplot(2,1,1);
 plot(v_ci:delta_v:v_co,lambda_stall);
 title('\lambda value at each wind speed');
 grid on;
 xlabel('Wind speed (m/s)');
 ylabel('\lambda');
 subplot(2,1,2);
 plot(v_ci:delta_v:v_co,cp_stall);
 title('C_{p} value at each wind speed');
 grid on;
 xlabel('Wind speed (m/s)');
 ylabel('C_{p}');
 
 %% Constant speed - Pitch regulated
 
 delta_v = 1;
 N = (v_co-v_ci)/delta_v + 1;
 
 speed = v_ci;
 
 N_beta = 26; % number of pitch angles that we will be evaluated
 col_betai = 13; % column of the first value of beta evaluated
 
 lambda_cspr = lambda_stall;
 cp_cspr = zeros(N_beta,N);
 cp_pr = zeros(1,N);
 Power_cspr = zeros(1,N);
 beta_pr = zeros(1,N);
 Power_cspr_rated = zeros(1,N); 
 cp_pcrated = zeros(1,N);
 Power_pc_rated = zeros(1,N);
 for i = 1:N
     
     for j = 1:N_beta
         n_values = 23;
         
         for ii = 1:23
             if isnan(performance(24-ii,col_betai + 4*(j-1)))
                 n_values = n_values - 1;
             end
         end
         
         col_cp = col_betai + 2 + 4*(j-1);
         col_lambda = col_cp -1;
         
         cp_cspr(j,i) = interp1(performance(1:n_values,col_lambda),performance(1:n_values,col_cp),lambda_stall(i));
         
         if cp_cspr(j,i) >= cp_pr(i)
             cp_pr(i) = cp_cspr(j,i);
             beta_pr(i) = performance(1,col_betai + 4*(j-1));
         end
     end
     
     Power_cspr(i) = 1 / 2 * cp_pr(i) * rho_air * pi * R_wt^2 * speed ^3;
     Power_cspr_rated(i) = Power_cspr(i);
     
     if Power_cspr_rated(i) > PR_wt * 10^6
         Power_cspr_rated(i) = PR_wt * 10^6;
         cp_localmax = 2 * PR_wt * 10^6 / (rho_air * speed ^3 * pi * R_wt^2);
         
         for j=1:N_beta
             if cp_cspr(j,i) <= cp_localmax
                 if cp_cspr(j,i) > cp_pcrated(i)
                     cp_pcrated(i) = cp_cspr(j,i);
                 end
             end
         end
         
     else
         cp_pcrated(i) = cp_pr(i);
     end
     
     Power_pc_rated(i) = cp_pcrated(i) / 2 * rho_air * pi * R_wt^2 * speed ^3;
     
    speed = speed + delta_v; 
    
 end
 
 figure; 
 plot(v_ci:delta_v:v_co,Power_cspr_rated,'DisplayName','Rated power');
 title('Power output of constant speed - pitch regulated wind turbine')
 grid on;
 xlabel('Wind speed (m/s)');
 ylabel('Power (W)');
 hold on;
 plot(v_ci:delta_v:v_co,Power_pc_rated,'DisplayName','Constant speed - pitch regulated');
 legend;
 
 
 figure;
 subplot(2,1,1);
 plot(v_ci:delta_v:v_co,lambda_stall);
 title('\lambda value at each wind speed');
 grid on;
 xlabel('Wind speed (m/s)');
 ylabel('\lambda');
 subplot(2,1,2);
 plot(v_ci:delta_v:v_co,cp_pr);
 title('C_{p} value at each wind speed');
 grid on;
 xlabel('Wind speed (m/s)');
 ylabel('C_{p}');
 
 figure; 
 plot(v_ci:delta_v:v_co,beta_pr);
 title('Pitch angle variation for a constant speed - pitch regulated wind turbine');
 grid on;
 xlabel('Wind speed (m/s)');
 ylabel('\beta (\textordmasculine)');

 %% VARIABLE SPEED - PITCHCONTROL
 
omega_max = V_tip / R_wt;
 
v_omegamax = V_tip / lambda_cpmax;
 
speed = v_ci;

cp_opt = zeros(1,N);
lambda_max = zeros(1,N);
beta_opt = zeros(1,N);
lambda_opt = zeros(1,N);
power_pc = zeros(1,N);
omega = zeros(1,N);
cp_prvs = zeros(N_beta,N);
cp_vspr = zeros(1,N);
Power_vspr = zeros(1,N);
Power_vspr_rated = zeros(1,N);
beta_vspr = zeros(1,N);

 for i = 1:N
     
     if speed <= v_omegamax
        Power_vspr(i) = 1 / 2 * rho_air * cp_max * pi * R_wt^2 * speed^3;
        Power_vspr_rated(i) = Power_vspr(i);
        lambda_opt(i) = lambda_cpmax;
        lambda_max(i) = lambda_cpmax;
        beta_opt(i) = beta_cpmax;
        cp_opt(i) = cp_max;
        omega(i) = lambda_opt(i) * speed / R_wt;
     elseif speed <= V_R
         
         lambda_max(i) = V_tip / speed;
         lambda_opt(i) = lambda_max(i);
         for j = 1:N_beta
             
         n_values = 23;
         
         for ii = 1:23
             if isnan(performance(24-ii,col_betai + 4*(j-1)))
                 n_values = n_values - 1;
             end
         end
         
         col_cp = col_betai + 2 + 4*(j-1);
         col_lambda = col_cp -1;
         
         cp_prvs(j,i) = interp1(performance(1:n_values,col_lambda),performance(1:n_values,col_cp),lambda_stall(i));
         
         if cp_prvs(j,i) >= cp_vspr(i)
             cp_vspr(i) = cp_prvs(j,i);
             beta_vspr(i) = performance(1,col_betai + 4*(j-1));
         end
         
         end
     
     Power_vspr(i) = 1 / 2 * cp_vspr(i) * rho_air * pi * R_wt^2 * speed ^3;
     Power_vspr_rated(i) = Power_vspr(i);
     omega(i) = lambda_opt(i) * speed / R_wt;
        
     else 
         
         lambda_max(i) = V_tip / speed;
         cp_local_max = PR_wt * 10^6 * 2 / (rho_air * pi * R_wt^2 * speed^3);
         
         for j = 15:4:115  %iteration columns
             for iii = 1:23 %iteration rows
                if isnan(performance(iii,j))
                    
                elseif performance(iii,j) < cp_local_max
                    if performance(iii,j-2) < beta_opt(i-1)
                        
                    elseif performance(iii,j) > cp_opt(i)
                        if performance(iii,j-1) < lambda_max(i)
                            cp_opt(i) = performance(iii,j);
                            beta_opt(i) = performance(1,j-2);
                            lambda_opt(i) = performance(iii,j-1);
                        end
                    end
                end
             end
         end
         omega(i) = lambda_opt(i) * speed / R_wt;
         Power_vspr(i) = 1 / 2 * rho_air * cp_opt(i) * pi * R_wt^2 * speed^3;
         Power_vspr_rated(i) = PR_wt * 10^6;
        
     end
     
     speed = speed + delta_v;

 end
 
 figure;
 plot(v_range2,beta_opt);
 grid on;
 title('Variable speed - Pitch control');
 xlabel('Wind speed (m/s)');
 ylabel('\beta (deg)');
 
 figure;
 plot(v_ci:delta_v:v_co,Power_vspr_rated,'DisplayName','Rated power');
%  line([v_omegamax,v_omegamax],[0,657065],'color','black','linestyle','--');
%  line([V_R,V_R],[0,1.43*10^6],'color','black','linestyle','--');
 title('Power output variable speed - pitch regulated wind turbine');
 xlabel('Wind speed (m/s)');
 ylabel('Power (W)');
 grid on;
 hold on;
 plot(v_ci:delta_v:v_co,Power_vspr,'DisplayName', 'Variable speed - Pitch regulated');
 plot(v_ci:delta_v:v_co,Power_pc_rated,'DisplayName','Constant speed - pitch regulated');
%  plot(v_ci:delta_v:v_co,Power_stall,'DisplayName','Stall wind turbine');
 % hold on;
 plot(v_range2,Power_stall_rated_2,'DisplayName', 'Constant speed - Stall regulated');
 legend;
 
 figure;
 subplot(2,1,1);
 plot(v_range2,lambda_opt,'DisplayName','\lambda wind turbine');
 hold on;
 plot(v_range2,lambda_max,'DisplayName','\lambda max');
 grid on;
 title('\lambda values');
 xlabel('Wind speed (m/s)');
 ylabel('\lambda');
 legend;
 
 subplot(2,1,2);
 plot(v_range2,omega);
 grid on;
 title('\Omega value for Variable Speed - Pitch Regulated wind turbine');
 xlabel('Wind speed (m/s)');
 ylabel('\Omega (rad/s)');

 %% ENERGY CALCULATION
 
 i = 0.06; % discount rate
 N = 20; % years
 R = i / (1-(1+i)^-N);
 
 % Average annual power output
 
    % Stall wind turbine
    sum_power = 0;
    N = (v_co-v_ci)/delta_v + 1;
    for i=2:N-1
        speed = v_ci + (i-1)*delta_v;
        if speed > V_R
        else
        sum_power = sum_power + weibull(speed, k, c) * Power_stall_rated_2(i);
        end
    end
    
    AAP_stall = delta_v * ((weibull(v_co, k, c)*0*Power_stall_rated_2(N) - weibull(v_ci, k, c)*Power_stall_rated_2(1))/2 +sum_power);
    AAP_stall_dens = AAP_stall / A_rotor;
    AEO_stall = AAP_stall * 8760;
    CF_stall = AAP_stall / (PR_wt*10^6);
    CAPEX_stall = (700 + 600) * PR_wt *10^(3);
    LCOE_stall = (CAPEX_stall * (R + 0.03))/ (AEO_stall * 10^-3);
    
    % Constant speed - pitch regulated
    sum_power = 0;
    sum_power2 = 0;
    N = (v_co-v_ci)/delta_v + 1;
    for i=2:N-1
        speed = v_ci + (i-1)*delta_v;
        sum_power = sum_power + weibull(speed, k, c)*Power_cspr_rated(i);
        sum_power2 = sum_power2 + weibull(speed, k, c)*Power_pc_rated(i);
    end
    
    AAP_pr = delta_v * ((weibull(speed, k, c)*Power_cspr_rated(N) -...
                   weibull(speed, k, c)*Power_cspr_rated(1))/2 +sum_power);
    AAP_pc = delta_v * ((weibull(speed, k, c)*Power_pc_rated(N) -...
                   weibull(speed, k, c)*Power_pc_rated(1))/2 +sum_power2);
    AAP_pr_dens = AAP_pr / A_rotor;
    AAP_pc_dens = AAP_pc / A_rotor;
    AEO_pr = AAP_pr * 8760;
    AEO_pc = AAP_pc * 8760;
    CF_pc = AAP_pc / (PR_wt*10^6);
    CF_pr = AAP_pr/ (PR_wt*10^6);
    CAPEX_pc = (900 + 600) * PR_wt *10^(3);
    LCOE_pc = (CAPEX_pc * (R + 0.03))/ (AEO_pr * 10^-3);
    
    % Variable speed - pitch regulated
    sum_power = 0;
    sum_power2 = 0;
    N = (v_co-v_ci)/delta_v + 1;
    for i=2:N-1
        speed = v_ci + (i-1)*delta_v;
        sum_power = sum_power + weibull(speed)*Power_vspr(i);
        sum_power2 = sum_power + weibull(speed)*Power_vspr_rated(i);
    end
    
    AAP_vspr = delta_v * ((weibull(speed)*Power_vspr(N) - weibull(speed)*Power_vspr(1))/2 +sum_power);
    AAP_vspr_rated = delta_v * ((weibull(speed)*Power_vspr_rated(N) - weibull(speed)*Power_vspr_rated(1))/2 +sum_power2);
    AAP_vspr_dens = AAP_vspr / A_rotor;
    AAP_vspr_rated_dens = AAP_vspr_rated / A_rotor;
    AEO_vspr = AAP_vspr_rated * 8760;
    CF_vspr = AAP_vspr_rated / (PR_wt * 10^6);
    CAPEX_vspr = (950 + 600)* PR_wt *10^(3);
    LCOE_vspr = (CAPEX_vspr * (R + 0.03))/ (AEO_vspr*10^-3);
