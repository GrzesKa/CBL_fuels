function plotKPIBarGraphs(injections, Efficiency, BSCO2, BSNOx, BSFC)
% Define the loads and efficiency values
loads = [2, 3, 4];        % Numeric values for loads (in bar)

for i=1:3
load = i+1;
% Creates the bar graph for the Efficiency KPI
figure;
subplot(2, 2, 1);
bar(injections, Efficiency(i,:), 'FaceColor', [0.2, 0.6, 0.8]);

% Add labels and title
%xlabel('Loads (bar)');
xlabel('Injection timing (-)');
ylabel('Efficiency (-)');

ylim([0,1]); % Adjust y-axis limits
title('Efficiency');

% Optional: Adjust appearance
grid on; % Add a grid for better readability
    ax = gca; % Get current axes
    ax.GridAlpha = 0.3; % Set grid line transparency
    ax.LineWidth = 1.2; % Make axis lines slightly thicker
    ax.FontSize = 10; % Adjust font size for axis labels and ticks
    ax.FontWeight = 'bold'; % Make axis labels bold


% Creates the bar graph of BSCO2 KPI
subplot(2, 2, 2);
bar(injections, BSCO2(i,:), 'FaceColor', [1, 0.5, 0]);

% Add labels and title
%xlabel('Loads (bar)');
xlabel('Injection timing (-)');
ylabel('BSCO2 (g/KWhr)');


title('BSCO2');

% Optional: Adjust appearance
grid on; % Add a grid for better readability
    ax = gca; % Get current axes
    ax.GridAlpha = 0.3; % Set grid line transparency
    ax.LineWidth = 1.2; % Make axis lines slightly thicker
    ax.FontSize = 10; % Adjust font size for axis labels and ticks
    ax.FontWeight = 'bold'; % Make axis labels bold


% Creates the bar graph of BSNOx KPI
subplot(2, 2, 3);
bar(injections, BSNOx(i,:), 'FaceColor', [0.4, 0.7, 0.3]);

% Add labels and title
%xlabel('Loads (bar)');
xlabel('Injection timing (-)');
ylabel('BSNOx (mg/KWhr)');

%ylim([0,1]); % Adjust y-axis limits
title('BSNox');

% Optional: Adjust appearance
grid on; % Add a grid for better readability
    ax = gca; % Get current axes
    ax.GridAlpha = 0.3; % Set grid line transparency
    ax.LineWidth = 1.2; % Make axis lines slightly thicker
    ax.FontSize = 10; % Adjust font size for axis labels and ticks
    ax.FontWeight = 'bold'; % Make axis labels bold



% Creates the bar graph of BSFC KPI
subplot(2, 2, 4);
bar(injections, BSFC(i,:), 'FaceColor', [1, 0, 0]);

% Add labels and title
%xlabel('Loads (bar)');
xlabel('Injection timing (-)');
ylabel('BSCFC (g/KWhr)');

%ylim([0,1]); % Adjust y-axis limits
title('BSFC');

% Optional: Adjust appearance
grid on; % Add a grid for better readability
    ax = gca; % Get current axes
    ax.GridAlpha = 0.3; % Set grid line transparency
    ax.LineWidth = 1.2; % Make axis lines slightly thicker
    ax.FontSize = 10; % Adjust font size for axis labels and ticks
    ax.FontWeight = 'bold'; % Make axis labels bold
sgtitle(sprintf('KPIs for load %d', load), 'FontSize', 12, 'FontWeight', 'bold');
hold off;
end