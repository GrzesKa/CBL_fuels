function plotAveragePressure(V_matrix, smooth_p, dm, bara, iselect)
    % PLOTSELECTEDCYCLE - Plots the p-V diagram for a selected cycle in both linear and log-log scales.
    %
    % Inputs:
    % V_matrix - Volume matrix [m^3]
    % p_matrix - Pressure matrix [Pa]
    % dm - Conversion factor to dm^3
    % bara - Conversion factor to bar
    % iselect - Index of the selected cycle

    % Extract data for the selected cycle
    V_cycle = V_matrix(:, iselect);
    p_cycle = smooth_p

    % Create and configure figure
    f2 = figure;
    set(f2, 'Position', [200, 400, 600, 800]); % Custom figure size

    % Subplot 1: Linear scale
    subplot(2, 1, 1);
    plot(V_cycle / dm^3, p_cycle / bara, 'LineWidth', 1);
    xlabel('Volume [dm^3]');
    ylabel('Pressure [bar]');
    xlim([0 0.8]);
    ylim([0.5 50]);
    set(gca, 'XTick', 0:0.1:0.8, 'XGrid', 'on', 'YGrid', 'on');
    title({'p-V Diagram (Linear Scale)'});

    % Subplot 2: Log-log scale
    subplot(2, 1, 2);
    loglog(V_cycle / dm^3, p_cycle / bara, 'LineWidth', 1);
    xlabel('Volume [dm^3]');
    ylabel('Pressure [bar]');
    xlim([0.02 0.8]);
    ylim([0.5 50]);
    set(gca, 'XTick', [0.02 0.05 0.1 0.2 0.5 0.8], ...
             'YTick', [0.5 1 2 5 10 20 50], ...
             'XGrid', 'on', 'YGrid', 'on');
    title({'p-V Diagram (Log-Log Scale)'});
end
