function [Efficiency_all, BSCO2_all, BSNOx_all, BSFC_all] = KPI_function(V_cycle, W_per_cycle,CO2, NOx,VolumeEmission)

RPM = 1500 ; % Rotations Per Minute of engine
Power_engine = ((W_per_cycle/1000) * RPM)/120 ; % Power of engine, work hardcoded due to errors
LHV_B7 = 43e6; % LHV of B7 diesel in J/kg
V_air =  max(V_cycle); % volume of air in 0.5 cycle * nr of cycles in 1 second
mass_fuel = 0.0001176; % mass of fuel being injected (data from sensors)
ro_air = 1.293 ; % Density air 
V_air_ps = V_air*RPM/120 ; % Volume of air per second
m_air = V_air_ps * ro_air ; % mass of air entering, calculated using v of air in 1 second * density
exhaust_massflow = mass_fuel + m_air; % total mass flowing
exhaust_massflow_grams = exhaust_massflow * 1000; % used for the brake specific KPI's in g/
Efuel = mass_fuel*LHV_B7 ; 

%% KPIs calculation 

Efficiency_all = Power_engine*1000/Efuel ; % work per cycle/Efuel (work in J)

%Molar mass
NOx_NO2 = 0.1; % Percentage of NOx that transforms into NO2
molar_CO2 = 44;   % g/mol
molar_NOx = (30 * (1 - NOx_NO2) + 46 * NOx_NO2); % Weighted molar mass of NOx (g/mol)
molar_H2O = 18;   % g/mol
molar_Ar = 40;
molar_N2 = 28;    % g/mol
molar_O2 = 32;    % g/mol


% Combustion volume fractions of diesel
emission_N2 = 0.76;         % volume fraction for N2
emission_O2 = 0.135;        % volume fraction for O2
emission_CO2 = 0.0525;      % volume fraction for CO2
emission_H20 = 0.05;        % volume fraction for H2O
emission_Ar = 0.008;
emission_NOx = 0.00125;     % volume fraction for NOx


% Total molar mass contribution
total_molar_mass = (emission_N2 * molar_N2) + (emission_O2 * molar_O2) + (emission_CO2 * molar_CO2) + (emission_H20 * molar_H2O) + (emission_Ar * molar_Ar) + (emission_NOx * molar_NOx);


% Volume-to-mass fraction conversion
mass_fraction_CO2 = (emission_CO2 * molar_CO2) / total_molar_mass;
mass_fraction_NOx = (emission_NOx * molar_NOx) / total_molar_mass;

% Mass flow rates OLD
%CO2_massflow = exhaust_massflow_grams * mass_fraction_CO2; 
%NOx_massflow = exhaust_massflow_grams * mass_fraction_NOx; 

% Mass flow rates NEW
DensityCO2 = 1.98e3;                                %Density of CO2 at 20 degrees celcius (g/m^3)

CO2_massflow = CO2/100*VolumeEmission*DensityCO2/100;
NOx_massflow = NOx*46.01/24*VolumeEmission/100;

%KPIs - x3600 to get to kWhr
BSCO2_all = (CO2_massflow*3600)/Power_engine;
BSNOx_all = (NOx_massflow*3600)/Power_engine;
BSFC_all = (mass_fuel*3600)/Power_engine;

disp(['Efficiency_all:', num2str(Efficiency_all)]);

disp(['BSCO2_all:', num2str(BSCO2_all), 'g/kWhr']);

disp(['BSNOx_all:', num2str(BSNOx_all), 'g/kWhr']);

disp(['BSFC_all:', num2str(BSFC_all), 'g/kWhr']);


end