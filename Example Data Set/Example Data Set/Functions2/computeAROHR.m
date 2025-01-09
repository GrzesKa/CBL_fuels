function [aROHR_all, aROHR] = computeAROHR(Ca, Gamma_at_angle, smooth_p, V_avg, smooth_dVdCa, smooth_dpdCa)
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
gamma = 1.2;
aROHR_all = zeros(1, length(gamma_full)); % Preallocate for speed, assuming numDatasets is defined
for i = 1:length(gamma_full)
    gammaidk = gamma_full(i);
    smoothpidk = smooth_p(i);
    Vavgidk = V_avg(i);
    smoothdpdCaidk = smooth_dpdCa(i);
    smoothdVdCaidk = smooth_dVdCa(i);
    % Compute aROHR using the interpolated gamma and inputs
    aROHR_all(i) = (gammaidk ./ (gammaidk - 1)) .* smoothpidk .* smoothdVdCaidk + (1 ./ (gammaidk - 1)) .* Vavgidk .* smoothdpdCaidk;

end


gamma = 1.4;
aROHR = (gamma / (gamma - 1)) * smooth_p .* smooth_dVdCa + (1 / (gamma - 1)) * V_avg .* smooth_dpdCa;