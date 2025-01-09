function [W_per_cycle, m_per_cycle, V_avg, mass_flow_fuel_cycle, AFR] = calculateMassflow(Cyl, RPM, selectedFuel, FuelTable, V_cycle, m_fuel_cycle, W_all)
    % Function to calculate cycle metrics such as work per cycle, air-fuel ratio, etc.

    % Inputs:
    % Cyl - Cylinder dimensions (structure containing Bore and Stroke)
    % RPM - Engine revolutions per minute
    % selectedFuel - Selected fuel type
    % FuelTable - Table containing fuel properties (e.g., Density)
    % V_cycle - Volume per cycle
    % m_fuel_cycle - Fuel mass per cycle
    % W_all - Array of work per cycle
    % intake_pressure - Intake pressure in Pascals
    % intake_temp - Intake temperature in Kelvin
    % spec_gasct_air - Specific gas constant for air

    % Outputs:
    % W_per_cycle - Average work per cycle
    % m_per_cycle - Average mass per cycle
    % V_avg - Average volume per cycle
    % mass_flow_fuel_cycle - Mass flow of fuel per cycle
    % AFR - Air-fuel ratio
    intake_pressure = 100000; %Pa ambient pressure 1 bar
    intake_temp = 293.15; %Kelvin ambient temperature
    spec_gasct_air = 287; %J/mol*K

    % Displaced volume per cycle
    volume_displaced_cycle = (pi / 4) * (Cyl.Bore)^2 * Cyl.Stroke;

    % Volumetric flow rate
    volumetric_flow = volume_displaced_cycle * 0.5 * (RPM / 60);

    % Air density
    density_air = intake_pressure / (intake_temp * spec_gasct_air);

    % Mass flow rate of air
    mass_flow_air = volumetric_flow * density_air;

    % Fuel density
    Densityrow = strcmp(FuelTable.Fuel, selectedFuel);
    densityfuel = FuelTable.Density(Densityrow);

    % Air mass per cycle
    m_air = (V_cycle - (m_fuel_cycle / densityfuel)) * density_air;

    % Air-fuel ratio
    AFR = m_air ./ m_fuel_cycle;

    % Mass flow of fuel
    mass_flow_fuel = mass_flow_air ./ AFR;

    % Mass flow of fuel per cycle
    mass_flow_fuel_cycle = mass_flow_fuel / 12.5;

    % Average volume per cycle
    V_avg = sum(V_cycle) / numel(W_all);

    % Cumulative work and fuel mass
    W_cumm = sum(W_all);
    m_fuel_cumm = sum(m_fuel_cycle);

    % Average mass and work per cycle
    m_per_cycle = m_fuel_cumm / numel(W_all);
    W_per_cycle = W_cumm / numel(W_all);

    % Display or save results
    disp('Work per cycle:');
    disp(W_per_cycle);
    disp('m per cycle:');
    disp(m_per_cycle);
    disp('Volume avg:');
    disp(V_avg);
end
