function [V_cycle] = plotPVDiagram(V_matrix, p_matrix, cycleNumber)

    if cycleNumber < 1 || cycleNumber > size(V_matrix, 2)
        error('Invalid cycle number. Please provide a cycle number between 1 and %d.', size(V_matrix, 2));
    end

    % Extract data for the specified cycle
    V_cycle = V_matrix(:, cycleNumber);
    p_cycle = p_matrix(:, cycleNumber);

    % Plot the p-V diagram
    figure;
    plot(V_cycle, p_cycle, 'LineWidth', 1.5);
    xlabel('Volume [m^3]');
    ylabel('Pressure [Pa]');
    title(sprintf('p-V Diagram for Cycle %d', cycleNumber));
    grid on;

    % Improve aesthetics if needed
    set(gca, 'FontSize', 12);
end
