function [W_all, V_avg, W_cumm, W_per_cycle, m_per_cycle] = calculateWork(V_matrix, p_matrix, m_fuel_matrix, Ncycles)
    % Function to calculate work per cycle, cumulative work, average work,
    % and average mass per cycle.
    %
    % Inputs:
    %   V_matrix     - Matrix of volumes (columns represent cycles)
    %   p_matrix     - Matrix of pressures (columns represent cycles)
    %   m_fuel_matrix- Matrix of fuel masses (columns represent cycles)
    %   Ncycles      - Number of cycles
    %
    % Outputs:
    %   W_all        - Work for each cycle
    %   V_avg        - Average volume over all cycles
    %   W_cumm       - Cumulative work over all cycles
    %   W_per_cycle  - Average work per cycle
    %   m_per_cycle  - Average mass of fuel per cycle

W_all = zeros(1, Ncycles);              % Preallocate work matrix
for i = 1:Ncycles
    V_cycle = V_matrix(:,i);
    p_cycle = p_matrix(:,i);
    W_all(i) = trapz(V_cycle, p_cycle); % Numerical integration (trapezoidal) to find work per cycle
    m_fuel_cycle = m_fuel_matrix(:,i); % Mass fuel for cycle i
end
V_avg = sum(V_cycle)/Ncycles;
W_cumm = sum(W_all); %Cummulative function of the Work
m_fuel_cumm = sum(m_fuel_cycle);
m_per_cycle = m_fuel_cumm/Ncycles;
W_per_cycle = W_cumm/Ncycles; %Average work
end

