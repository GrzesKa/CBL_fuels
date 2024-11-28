function [W] = ComputeW(p, V)
    %Inputs:
    % p_cycle : pressure [Pa]
    % V_cycle : volume [m^3]
    % 
    % Output:
    % W : work [J]
    
    % Step 1: Calculate work from p-V Diagram
    W = trapz(V,p,1);

end