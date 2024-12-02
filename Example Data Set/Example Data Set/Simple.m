%% Info
% Toerental: 1500 RPM
% SOA van 4.2º voor TDC
% Resolutie van 0.2º CA
% Data voor 69 cycles (maximale van de Smetec, de OGO gensets kunnen in principe “onbeperkt” aan)
% 
%% init
clear all; clc;close all;
addpath( "Functions","Nasa");
%% Units
mm      = 1e-3;
dm      = 0.1;
bara    = 1e5;
MJ      = 1e6;
kWhr    = 1000*3600;
volperc = 0.01; % Emissions are in volume percentages
ppm     = 1e-6; % Some are in ppm (also a volume- not a mass-fraction)
g       = 1e-3;
s       = 1;

Ncycles = 69;  % Define the number of cycles
RPM = 1500; %Defining the RPM
%% Fuel properties

% Select fuel
selectedFuel = 'Diesel';


% Define fuel properties
FuelTable = table(...
    {'Diesel'; 'HVO'; 'FAME'; 'GTL'}, ...        % Fuel names
    [43e6; 43.7e6; 38.3e6; 44e6], ...            % Lower Heating Value [J/kg]
    [52.2; 74.5; 55.2; 74], ...                  % Cetane number [-]
    [2.62; 1.5; 2.1; 2.5], ...                   % CO2 [kg/L]
    [9; 7; 11; 7], ...                           % NOx [g/kWh]
    [836.1; 764; 882; 777.1], ...                % Density [kg/m^3]
    [2.7638; 2.88; 4.43; 2.5774], ...            % Viscosity [mm^2/s]
    'VariableNames', {'Fuel','LHV', 'Cetane', 'CO2','NOx', 'Density','Viscosity'});

%% Load NASA maybe you need it at some point?
% Global (for the Nasa database in case you wish to use it).

global Runiv
Runiv = 8.314;
[SpS,El] = myload('Nasa\NasaThermalDatabase.mat',{'Diesel','O2','N2','CO2','H2O'});

%% Engine geom data (check if these are correct)
Cyl.Bore                = 104*mm;
Cyl.Stroke              = 85*mm;
Cyl.CompressionRatio    = 21.5;
Cyl.ConRod              = 136.5*mm;
Cyl.TDCangle            = 180;

% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).

CaIVO = -355;  %Intake valve opens
CaIVC = -135;  %Intake valve closes
CaEVO = 149;   %Exhaust valve opens
CaEVC = -344;  %Exhaust valve closes
CaSOI = -3.2;  %Start of Injection - CHANGE IF WE PLAY AROUND WITH IT IN THE EXPERIMENT

% Write a function [V] = CylinderVolume(Ca,Cyl) that will give you Volume
% for the given Cyl geometry. If you can do that you can create pV-diagrams

%% Load data (of txt file)
FullName            = fullfile('Data','ExampleDataSet.txt');
dataIn              = table2array(readtable(FullName));
[Nrows,Ncols]       = size(dataIn);                    % Determine size of array
NdatapointsperCycle = 720/0.2;                         % Nrows is a multitude of NdatapointsperCycle
Ncycles             = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
% Check the size of matrix matches
disp(['Rows in dataIn: ', num2str(Nrows)]);
disp(['Expected rows (NdatapointsperCycle * Ncycles): ', num2str(NdatapointsperCycle * Ncycles)]);

Ca_matrix           = reshape(dataIn(:,1),NdatapointsperCycle,Ncycles);      % Both p and Ca are now matrices of size (NCa,Ncycles)
p_matrix            = reshape(dataIn(:,2),NdatapointsperCycle,Ncycles)*bara; % This forms the Matrix For the p values for each Ca per cycle and converts the pressure to SI units 
m_fuel_matrix       = reshape(dataIn(:,3),NdatapointsperCycle,Ncycles);      % Same thing as for pressure, but without the bara function

%% Work and Volume Calculation

% Volume calculation
V_matrix = CylinderVolume(Ca_matrix, Cyl);  % Calculate volume at each crank angle

% Work calculation
W_all = zeros(1, Ncycles);              % Preallocate work matrix
for i = 1:Ncycles
    V_cycle = V_matrix(:,i);
    p_cycle = p_matrix(:,i);
    W_all(i) = trapz(V_cycle, p_cycle); % Numerical integration (trapezoidal) to find work per cycle
end

% Display or save results
disp('Work for each cycle:');
disp(W_all);

% Plot p-V diagram for one cycle (e.g., cycle 1)
figure;
plot(V_matrix(:,1), p_matrix(:,1));
xlabel('Volume [m^3]');
ylabel('Pressure [Pa]');
title('p-V Diagram for Cycle 1');
grid on;


% %% KPI function implementations
% 
% % BSFC calculation
% BSFC_all = zeros(1, Ncycles);  % empty matrix to store values
% 
% for i = 1:Ncycles
%     % Extract the pressure and volume data for the current cycle (i)
%     V_cycle = V_matrix(:, i);          % Volume for cycle i
%     p_cycle = p_matrix(:, i);          % Pressure for cycle i
%     m_fuel_cycle = m_fuel_matrix(:,i); % Mass fuel for cycle i
% 
%     % Calculate BSFC using the ComputeBSFC function
%     BSFC_all(i) = ComputeBSFC(p_cycle, V_cycle, RPM, m_fuel_cycle); % Calculates the BSFC values per cycle and adds them to the empty matrix
% end
% 
% disp('BSFC for each cycle:');
% disp(BSFC_all);
% 
% % Efficiency calculation
% 
% Efficiency_all = zeros(1, Ncycles); 
% 
% for i = 1:Ncycles
%     V_cycle = V_matrix(:, i);           % Volume for cycle i
%     p_cycle = p_matrix(:, i);           % Pressure for cycle i
%     m_fuel_cycle = m_fuel_matrix(:,i);  % Mass fuel for cycle i
% 
% 
%     Efficiency_all(i) = efficiency(m_fuel_cycle, selectedFuel, FuelTable, V_cycle, p_cycle, RPM); %Calculates the efficiency per cycle and adds to the empty matrix 
% end
% 
% disp('Efficiency for each cycle:');
% disp(Efficiency_all);
% 
% %CO2 Emissions (incomplete)
% 
% bsCO2_all = zeros(1, Ncycles);
% 
% for i = 1:Ncycles
%     m_fuel_cycle = m_fuel_matrix(:,i); % Mass fuel for cycle i
%     bsCO2(i) = CO2Emissions(m_fuel_matrix, P, selectedFuel, Fueltable);
% end
% 
% disp(bsCO2_all);
%% Plotting 
f1=figure(1);
set(f1,'Position',[ 200 800 1200 400]);             % Just a size I like. Your choice
pp = plot(Ca_matrix,p_matrix/bara,'LineWidth',1);   % Plots the whole matrix
xlabel('Ca');ylabel('p [bar]');                     % Always add axis labels
xlim([-360 360]);ylim([0 50]);                      % Matter of taste
iselect = 10;                                       % Plot cycle 10 again in the same plot to emphasize it. Just to show how to access individual cycles.
line(Ca_matrix(:,iselect),p_matrix(:,iselect)/bara,'LineWidth',2,'Color','r');
YLIM = ylim;
% Add some extras to the plot
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b'); % Plot a vertical line at IVC. Just for reference not a particular reason.
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r'); % Plot a vertical line at EVO. Just for reference not a particular reason.
set(gca,'XTick',[-360:60:360],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('All cycles in one plot.')

%% pV-diagram

% Extract data for the selected cycle
V_cycle = V_matrix(:, iselect);
p_cycle = p_matrix(:, iselect);

f2 = figure(2);
set(f2,'Position',[ 200 400 600 800]);              % Just a size I like. Your choice

subplot(2,1,1)
plot(V_matrix/dm^3,p_matrix(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');                      % Always add axis labels
xlim([0 0.8]);ylim([0.5 50]);                              % Matter of taste
set(gca,'XTick',[0:0.1:0.8],'XGrid','on','YGrid','on');    % I like specific axis labels. Matter of taste
title({'pV-diagram'})

subplot(2,1,2)
loglog(V_matrix/dm^3,p_matrix(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0.02 0.8]);ylim([0 50]);                      % Matter of taste
set(gca,'XTick',[0.02 0.05 0.1 0.2 0.5 0.8],...
    'YTick',[0.5 1 2 5 10 20 50],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram'})

%% Constants
%Cp = CpNasa(T, )
%Cv = CvNasa(T, )

%% Compute Volume and Derivatives

% Define variables
Ca = Ca_matrix;       % Use all crank angle data
p = p_matrix;         % Use all pressure data
V = V_matrix;         % Use all volume data

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

%% Compute aROHR
gamma = 1.2; % Update gamma if needed

aROHR = (gamma / (gamma - 1)) * smooth_p .* smooth_dVdCa + (1 / (gamma - 1)) * V_avg .* smooth_dpdCa;

%% Plot aROHR
f3 = figure(3);
set(f3, 'Position', [400 400 800 400]); % Figure size
plot(Ca(:, 1), aROHR, 'LineWidth', 1); % Plot averaged aROHR
xlabel('Crank Angle (°)');
ylabel('aROHR [J/°CA]');
xlim([-45 135]);
grid on;
title('Apparent Rate of Heat Release (aROHR) vs Crank Angle');

%% Find the Indices for CaSOI and CaEVO
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

%% Perform Cumulative Integration
aHR = cumtrapz(Ca_from_start, aROHR_from_start); % Cumulative heat release

%% Plot Cumulative Heat Release (aHR)
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
