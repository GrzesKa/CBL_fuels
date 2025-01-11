function plotAHR(Ca_from_start, aHR, CaSOI, idx_BDC, p)
figure;
plot(Ca_from_start, aHR, 'b', 'LineWidth', 2); % Plot cumulative heat release
hold on;

% Annotate 10%, 50%, and 90% heat release points
[~, idx10] = min(abs(aHR - 0.1 * max(aHR)));
[~, idx50] = min(abs(aHR - 0.5 * max(aHR)));
[~, idx90] = min(abs(aHR - 0.9 * max(aHR)));

plot(Ca_from_start(idx10), aHR(idx10), 'yo', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow'); % 10%
text(Ca_from_start(idx10), aHR(idx10) + 20, '10%', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

plot(Ca_from_start(idx50), aHR(idx50), 'yo', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow'); % 50%
text(Ca_from_start(idx50), aHR(idx50) + 20, '50%', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

plot(Ca_from_start(idx90), aHR(idx90), 'yo', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow'); % 90%
text(Ca_from_start(idx90), aHR(idx90) + 20, '90%', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

% Add reference lines
line([0 0], ylim, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5); % Vertical line at 0° CA

% Add plot details
xlabel('Crank Angle [deg]');
ylabel('aHR [J]');
title('Cumulative Heat Release (aHR) vs Crank Angle');
xlim([CaSOI - 10, 130]); % Adjust x-axis limits
ylim([min(aHR) - 100, max(aHR) + 100]); % Adjust y-axis limits
grid on;
hold off;

pbdx = p(idx_BDC);