function [smooth_dVdCa, smooth_dpdCa, smooth_p, V_avg, Ca, p, V, Ca_single, NCa] = calculatingDerivatives(Ca_matrix, p_matrix, V_matrix, Ncycles)
% Define variables
Ca = Ca_matrix;       % Use all crank angle data
p = p_matrix;         % Use all pressure data
V = V_matrix;         % Use all volume data
Ca_single = Ca(:, 1); % Use the crank angle array (same for all cycles)
% Preallocate derivative matrices
dVdCa = zeros(size(V));
dpdCa = zeros(size(p));

% Compute derivatives for each cycle
for i = 1:Ncycles
    dVdCa(:, i) = gradient(V(:, i), Ca(:, i));
    dpdCa(:, i) = gradient(p(:, i), Ca(:, i));
end

% Average over all cycles and apply filtering
NCa = size(dVdCa, 1);  % number of unique crank angles

% Initialize averaged variables
filtered_averaged_dVdCa = zeros(NCa, 1); 
filtered_averaged_dpdCa = zeros(NCa, 1); 
filtered_averaged_p = zeros(NCa, 1);

for i = 1:NCa
    % dVdCa
    angle_dVdCa = dVdCa(i, :);
    sorted_dVdCa = sort(angle_dVdCa);
    n = length(sorted_dVdCa);
    lower_bound = ceil(n * 0.025);  
    upper_bound = floor(n * 0.975); 
    filtered_dVdCa = sorted_dVdCa(lower_bound:upper_bound);
    filtered_averaged_dVdCa(i) = mean(filtered_dVdCa);
    
    % dpdCa
    angle_dpdCa = dpdCa(i, :);
    sorted_dpdCa = sort(angle_dpdCa);
    filtered_dpdCa = sorted_dpdCa(lower_bound:upper_bound);
    filtered_averaged_dpdCa(i) = mean(filtered_dpdCa);
    
    % Pressure
    angle_p = p(i, :);
    sorted_p = sort(angle_p);
    filtered_p = sorted_p(lower_bound:upper_bound);
    filtered_averaged_p(i) = mean(filtered_p);
end

% Smooth the filtered data
smooth_dVdCa = sgolayfilt(filtered_averaged_dVdCa, 1, 9);
smooth_dpdCa = sgolayfilt(filtered_averaged_dpdCa, 1, 9);
smooth_p = sgolayfilt(filtered_averaged_p, 1, 9);

% Use the average volume (since geometry doesn't change)
V_avg = mean(V, 2); % Average over all cycles
