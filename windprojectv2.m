%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% WIND ANALYSIS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%  ________________________________________________________________________
%% DATA AND PARAMETERS 
% *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *

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
