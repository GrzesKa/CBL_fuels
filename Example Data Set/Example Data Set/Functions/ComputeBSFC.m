function [BSFC] = ComputeBSFC(p_cycle, V_cycle, RPM, m_fuel_cycle)
    % Inputs:
    % p      : Pressure data [Pa] (vector for one cycle)
    % V      : Volume data [m^3] (vector for one cycle)
    % RPM    : Engine speed [Revolutions Per Minute]
    % m_fuel  : Fuel mass flow rate [kg/s]
    % Cyl    : Cylinder geometry struct (for additional parameters if needed)
    %
    % Output:
    % BSFC   : Brake Specific Fuel Consumption [kg/kWh]
    
    % Step 1: Calculate Work from p-V Diagram
    Work = trapz(V_cycle, p_cycle); % Numerical integration of p-V curve (Work in J)
    
    % Step 2: Compute Power
    %N_revs_per_cycle = 2; % For a 4-stroke engine
    %Power = (Work * RPM) / (N_revs_per_cycle * 60); % Power in Watts (J/s)
   

    % Step 3: Calculate BSFC
    %BSFC = (m_fuel_cycle / Power) * 3600; % Convert to kg/kWh

    power_output = Work*RPM / 60;

    fuel_consumed = sum(m_fuel_cycle);

    BSFC = fuel_consumed / power_output;
end