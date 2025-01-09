function aROHR = computeAROHR(Ca, Gamma_at_angle, smooth_p, V_cycle, dVdCa, dpdCa)
% Interpolate gamma for the given crank angles in Ca(:, 1)

    gamma_start = Gamma_at_angle(1, 2);  % First gamma value
    gamma_end = Gamma_at_angle(end, 2); % Last gamma value

    gamma_full = zeros(size(Ca(:, 1)));  % Initialize the output array

    % Crank angles less than the first
    gamma_full(Ca(:, 1) < Gamma_at_angle(1, 1)) = gamma_start;

    % Crank angles greater than the last
    gamma_full(Ca(:, 1) > Gamma_at_angle(end, 1)) = gamma_end;

    % Crank angles within the range
    in_range = (Ca(:, 1) >= Gamma_at_angle(1, 1)) & (Ca(:, 1) <= Gamma_at_angle(end, 1));
    gamma_full(in_range) = interp1(Gamma_at_angle(:, 1), Gamma_at_angle(:, 2), Ca(in_range, 1), 'linear');

% % Initialize aROHR_all

aROHR = (gamma_full ./ (gamma_full - 1)) .* smooth_p .* dVdCa(:,1) + ...
        (1 ./ (gamma_full - 1)) .* V_cycle .*dpdCa;