function [smooth_P,dpdCa] = Pegging_dpdCa(smooth_p,NCa, Ca)
%% Pegging

desired_value = 101325;
% Calculate the difference between the original first value and the adjusted BDC value
idx_BDC = NCa/4;
diff = smooth_p(idx_BDC) - desired_value;
    
% Pegging at 1 bar at BDC
for i = 1:length(smooth_p)

    idx_BDC = NCa/4;
    diff = smooth_p(idx_BDC) - desired_value;
    
    % Adjust the current element by adding the calculated difference
    smooth_P(i) = smooth_p(i) - diff;
end
smooth_P=smooth_P(:);
Ca_single = Ca(:, 1);
dpdCa = zeros(size(Ca_single));
dpdCa = gradient(smooth_P, Ca_single);
dpdCa = sgolayfilt(dpdCa, 5, 37);
end