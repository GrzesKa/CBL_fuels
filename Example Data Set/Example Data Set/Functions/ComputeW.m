function [W] = ComputeW(p_cycle, V_cycle)
    %Inputs:
    % p_cycle : pressure [Pa]
    % V_cycle : volume [m^3]
    % 
    % Output:
    % W : work [J]
    
    % Step 1: Calculate work from p-V Diagram
    W = trapz(V_cycle,p_cycle,1);

end