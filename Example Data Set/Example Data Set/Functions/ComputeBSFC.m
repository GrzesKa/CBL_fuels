function [BSFC] = ComputeBSFC(p, V, RPM, m_dot)
    % Inputs:
    % p      : Pressure data [Pa] (vector for one cycle)
    % V      : Volume data [m^3] (vector for one cycle)
    % RPM    : Engine speed [Revolutions Per Minute]
    % m_dot  : Fuel mass flow rate [kg/s]
    % Cyl    : Cylinder geometry struct (for additional parameters if needed)
    %
    % Output:
    % BSFC   : Brake Specific Fuel Consumption [kg/kWh]
    
    % Step 1: Calculate Work from p-V Diagram
    Work = trapz(V, p); % Numerical integration of p-V curve (Work in J)
    
    % Step 2: Compute Power
    N_revs_per_cycle = 2; % For a 4-stroke engine
    Power = (Work * RPM) / (N_revs_per_cycle * 60); % Power in Watts (J/s)
    
    % Step 3: Calculate BSFC
    BSFC = (m_dot / Power) * 3600; % Convert to kg/kWh
end