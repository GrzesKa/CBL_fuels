function [Efficiency_all, BSCO2_all, BSNOx_all, BSFC_all] = KPI_function(V_cycle, W_per_cycle,CO2, NOx,VolumeEmission,FuelTable,selectedFuel,smooth_P)

RPM = 1500 ; % Rotations Per Minute of engine
massflow_fuel = 0.0002; % kg/s
kWhr_conversion = 3.6E6;
Power_engine = (W_per_cycle/kWhr_conversion);


%% KPIs calculation 

Efficiency_all = efficiency(massflow_fuel, selectedFuel, FuelTable, V_cycle, smooth_P, RPM);



% Mass flow rates NEW
DensityCO2 = 1.98e3;                                %Density of CO2 at 20 degrees celcius (g/m^3)

CO2_massflow = CO2/100*VolumeEmission*DensityCO2/100;
NOx_massflow = NOx*(46.01*0.15+30.1*0.85)/24*VolumeEmission/100;

%KPIs - x3600 to get to kWhr
BSCO2_all =CO2_massflow/Power_engine;
BSNOx_all = NOx_massflow/Power_engine;
BSFC_all = massflow_fuel/Power_engine;

disp(['Efficiency_all:', num2str(Efficiency_all)]);
disp(['BSCO2_all:', num2str(BSCO2_all), 'g/kWhr']);
disp(['BSNOx_all:', num2str(BSNOx_all), 'g/kWhr']);
disp(['BSFC_all:', num2str(BSFC_all), 'g/kWhr']);


end