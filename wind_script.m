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
R = 287.05;
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
    wind_MAST(:,1), wind_MAST(:,2), wind_MAST(:,3), wind_MAST(:,4),...
    wind_MAST(:,9), wind_MAST(:,10)];
wind_MAST2 = wind_MAST; NM = size(wind_MAST2, 1); NC = size(wind_MAST2, 2);
wind_MAST_aux = zeros(NM,NC);
% Correct the data from the meteorological mast 
counter = 1;
initial_date = wind_MAST(1, 3:5);
last_date    = wind_MAST(end, 3:5);

year_diff    = last_date(1) - initial_date(1);
day_diff     = last_date(3) - initial_date(3);
month_diff   = last_date(2) - initial_date(2);
if month_diff < 0; year_diff = year_diff - 1; end

limit_date = initial_date; limit_date(1) = initial_date(1) + year_diff;
limit_date_pt = limit_date; limit_date_pt(1) = limit_date_pt(1) - 1;

pressure_aux = zeros(NM,1); temperature_aux = zeros(NM,1);
density_aux = zeros(NM,1);

% First correct for the pressure and temperature
for ii = 1:NM
    if wind_MAST2(ii,3:5) == limit_date_pt
        break
    elseif wind_MAST2(ii,8) >= 0 
        pressure_aux(counter)    = wind_MAST2(ii,8);
        temperature_aux(counter) = wind_MAST2(ii,7);
        density_aux(counter) = pressure_aux(counter) * 100 / ...
            (R * (temperature_aux(counter) + 273.15));
        counter = counter + 1;
    else
        pressure_aux(counter) = pressure_aux(counter - 1);
        temperature_aux(counter) = temperature_aux(counter - 1);
        density_aux(counter) = pressure_aux(counter) * 100 / ...
            (R * (temperature_aux(counter) + 273.15));
        counter = counter + 1;
    end
end

pressure = pressure_aux(1:(counter-1)); 
temperature = temperature_aux(1:(counter-1));
density = density_aux(1:(counter-1));

mean_pressure = mean(pressure);
mean_temperature = mean(temperature);
mean_density = mean(density);

% Then correct for the other values
counter = 1;
for ii = 1:NM
    if wind_MAST2(ii, 3:5) == limit_date
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