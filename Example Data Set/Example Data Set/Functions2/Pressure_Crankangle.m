function Pressure_Crankangle(Ca_matrix, p_matrix, bara, CaIVC, CaEVO, iselect)
    % PLOTALLCYCLES - Plots pressure vs crank angle for all cycles and highlights one cycle.
    %
    % Inputs:
    % Ca_matrix - Crank angle matrix [degrees]
    % p_matrix - Pressure matrix [Pa]
    % bara - Conversion factor to bar
    % CaIVC - Crank angle at intake valve closing [degrees]
    % CaEVO - Crank angle at exhaust valve opening [degrees]
    % iselect - Index of the cycle to emphasize

    % Create and configure figure
    f1 = figure;
    set(f1, 'Position', [200, 800, 1200, 400]); % Custom figure size

    % Plot all cycles
    plot(Ca_matrix, p_matrix / bara, 'LineWidth', 1);
    hold on;

    % Highlight selected cycle
    plot(Ca_matrix(:, iselect), p_matrix(:, iselect) / bara, 'LineWidth', 2, 'Color', 'r');

    % Add axis labels and title
    xlabel('Crank Angle [°]');
    ylabel('Pressure [bar]');
    title('All Cycles in One Plot');

    % Set axis limits
    xlim([-360 360]);
    ylim([0 50]);

    % Add vertical lines for IVC and EVO
    YLIM = ylim;
    line([CaIVC, CaIVC], YLIM, 'LineWidth', 1, 'Color', 'b');
    line([CaEVO, CaEVO], YLIM, 'LineWidth', 1, 'Color', 'r');

    % Configure grid and ticks
    set(gca, 'XTick', -360:60:360, 'XGrid', 'on', 'YGrid', 'on');

    % Finalize
    hold off;
end
