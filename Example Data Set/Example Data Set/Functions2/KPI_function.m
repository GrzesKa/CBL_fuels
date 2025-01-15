function [Efficiency_all, BSCO2_all, BSNOx_all, BSFC_all] = KPI_function(V_cycle, W_per_cycle,CO2, NOx,VolumeEmission,FuelTable,selectedFuel,smooth_P)

RPM = 1500 ; % Rotations Per Minute of engine
massflow_fuel = 0.0002; % kg/s
kWhr_conversion = 3.6E6;
Power_engine = (W_per_cycle/kWhr_conversion);
s = 120 / RPM;

%% KPIs calculation 

Efficiency_all = efficiency(massflow_fuel, selectedFuel, FuelTable, V_cycle, smooth_P, RPM);



% Mass flow rates NEW
DensityCO2 = 1.98e3;                                %Density of CO2 at 20 degrees celcius (g/m^3)

CO2_massflow = CO2/100*VolumeEmission*DensityCO2/1000;
NOx_massflow = NOx*(46.01*0.15+30.1*0.85)/24*VolumeEmission/100;


BSCO2_all =CO2_massflow/Power_engine;
BSNOx_all = NOx_massflow/Power_engine;
BSFC_all = massflow_fuel*s*1000/Power_engine;


end