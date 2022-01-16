function output = wind_analysis(wind_data, figure_names, nbins)

    %% WIND ANALYSIS %%
    
    sz_wind = size(wind_data);
    N = sz_wind(1);
    V_mean = 0;
    v = zeros(N,1);
    
    %%%% Mean velocity
    for ii = 1:N
        v(ii) = wind_data(ii,1);
        V_mean = V_mean + v(ii);
    end
    
    V_mean = V_mean / N;
    
    %%%% Standard deviation
    sqdif = 0;
    for ii=1:N
      sqdif = sqdif + (v(ii) - V_mean)^2;
    end
    
    sigma = sqrt(sqdif / (N - 1)); %% standard deviation
    
    I = sigma / V_mean; %% turbulence intensity
    
    %% WIND ROSE
    
    d = zeros(sz_wind(1),1); %% Vector of the wind directions
    
    for ii=1:sz_wind(1)
        d(ii) = wind_data(ii,2);
    end
        
    
    Properties = {'anglenorth',0,'angleeast',90,...
        'labels',{'N (0)','S (180)','E (90)','W (270)'},...
        'freqlabelangle',22.5};
    WindRose(d,v,Properties);
%     set(gcf, 'PaperPosition', [0 0 15 10]); 
%     set(gcf, 'PaperSize', [15 10]); 
    print(figure_names.windrose,'-dpdf');

    %% WEIBULL DISTRIBUTION (semiempirical formulas of Bowden)
    
    c = 2 * V_mean / sqrt(pi); %scale parameter
    k = (sigma / V_mean)^(-1.086); %shape parameter
    
    figure;
    v_range = 0:0.01:30;
    p_v = k / c .* (v_range / c) .^(k-1) .* exp(-(v_range/c).^k);
    
    figure('Name', 'Weibull plot');
    plot(v_range,p_v, 'k-');
    xlabel('Wind speed [m/s]', 'Interpreter', 'Latex');
    ylabel('Probability f(v)', 'Interpreter', 'Latex');
    title('\textbf{Weibull probability distribution}','Interpreter', 'Latex');
    grid on;
    text(15,0.07,{['Scaling factor, A = ',num2str(c)],...
        ['Shape factor, k = ', num2str(k)]}, 'Interpreter', 'Latex')
    hold on; gcf;
    histogram(wind_data(:,1), nbins, 'Normalization','probability');
    set(gcf, 'PaperPosition', [0 0 15 10]); 
    set(gcf, 'PaperSize', [15 10]); 
    print(figure_names.weibull,'-dpdf');

    output.I = I; % Turbulence intensity

end