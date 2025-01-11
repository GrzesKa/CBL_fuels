function [Ca_from_start, aROHR_from_start] = indices(Ca, CaSOI, CaEVO, aROHR)

Ca_single = Ca(:, 1); % Use the crank angle array (same for all cycles)
[~, idx_start] = min(abs(Ca_single - CaSOI)); % Index closest to CaSOI
[~, idx_end] = min(abs(Ca_single - CaEVO));   % Index closest to CaEVO

% Ensure that idx_end is after idx_start
if idx_end <= idx_start
    error('CaEVO must be after CaSOI in the data');
end

% Slice data from CaSOI to CaEVO
Ca_from_start = Ca_single(idx_start:idx_end); % Crank angle from CaSOI to CaEVO
aROHR_from_start = aROHR(idx_start:idx_end);  % aROHR from CaSOI to CaEVO