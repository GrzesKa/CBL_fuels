%clear all; clc;close all;
addpath( "Functions","Nasa","Functions2");

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



%% Preallocate results structure
numLoads = 3;
%Experiment = 1
if Experiment == 1
    fuel = 'Diesel100';
    selectedFuel = 'Diesel'
elseif Experiment == 2
    fuel = 'GTL50';
    selectedFuel = 'GTL50'
elseif Experiment == 3
    fuel = 'GTL100';
    selectedFuel = 'GTL'
elseif Experiment == 4
    fuel = 'HVO50';
    selectedFuel = 'HVO50'
elseif Experiment == 5
    fuel = 'HVO100'
    selectedFuel = 'HVO'
end
DataSetTimings = DataSetCheck(fuel)

results = struct ( 'Work', [], 'BSFC', [], 'Efficiency', [], 'aROHR', [], 'CumulativeHR', []) ;

%% Fuel properties

% Select fuel
selectedFuel = 'Diesel';


% Define fuel properties
FuelTable = table(...
    {'Diesel'; 'HVO'; 'FAME'; 'GTL'; 'HVO50'; 'GTL50'}, ...                            % Fuel names
    [43e6; 43.7e6; 38.3e6; 44e6; mean([43e6, 43.7e6]); mean([43e6, 44e6])], ...        % Lower Heating Value [J/kg]
    [52.2; 74.5; 55.2; 74; mean([52.2, 74.5]); mean([52.2, 74])], ...                  % Cetane number [-]
    [2.62; 1.5; 2.1; 2.5; mean([2.62, 1.5]); mean([2.62, 2.5])], ...                   % CO2 [kg/L]
    [9; 7; 11; 7; mean([9, 7]); mean([9, 7])], ...                                     % NOx [g/kWh]
    [836.1; 764; 882; 777.1; mean([836.1, 764]); mean([836.1, 777.1])], ...            % Density [kg/m^3]
    [2.7638; 2.88; 4.43; 2.5774; mean([2.7638, 2.88]); mean([2.7638, 2.5774])], ...    % Viscosity [mm^2/s]
    'VariableNames', {'Fuel','LHV', 'Cetane', 'CO2','NOx', 'Density','Viscosity'});


%% Emission data

FullName = fullfile(sprintf('Data/FinalDataSet/%s data',fuel), sprintf('EmissionData%s.txt',fuel)); % Adjust file naming as needed
dataEmission = readtable(FullName);
dataEmission.Properties.VariableNames = {'Load', 'InjectionTiming', 'CO', 'CO2', 'HC', 'O2', 'NOx'};

%% Loop through all the data

for T = 1:length(DataSetTimings)
    timing = DataSetTimings(T);
% Emmision dataset sorting

% for (go through folder)
% Load dataset
%FullName            = fullfile('Data','ExampleDataSet.txt');
FullName = fullfile(sprintf('Data/FinalDataSet/%s data/%sPressureData',fuel,fuel), sprintf('%s_CA%d.txt',fuel,timing)); % choose one  file in folder, after executed choose next one etc.
dataIn = table2array(readtable(FullName));

[Nrows,Ncols]       = size(dataIn);                    % Determine size of array
NdatapointsperCycle = 720/0.2;                         % Nrows is a multitude of NdatapointsperCycle
Ncycles             = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on

Ca              = reshape(dataIn(:,1),[],Ncycles);
[V, Vmin, Vdiff] = CylinderVolume(Ca,Cyl);

% Check the size of matrix matches
% disp(['Rows in dataIn: ', num2str(Nrows)]);
% disp(['Expected rows (NdatapointsperCycle * Ncycles): ', num2str(NdatapointsperCycle * Ncycles)]);

Ca_matrix           = reshape(dataIn(:,1),NdatapointsperCycle,Ncycles);      % Both p and Ca are now matrices of size (NCa,Ncycles)
p_matrix            = reshape(dataIn(:,2),NdatapointsperCycle,Ncycles)*bara; % This forms the Matrix For the p values for each Ca per cycle and converts the pressure to SI units 
m_fuel_matrix       = reshape(dataIn(:,3),NdatapointsperCycle,Ncycles);      % Same thing as for pressure, but without the bara function


% Work and Volume Calculation

% Volume calculation
V_matrix = CylinderVolume(Ca_matrix, Cyl);  % Calculate volume at each crank angle

% Work calculation
W_all = zeros(1, Ncycles);              % Preallocate work matrix
for i = 1:Ncycles
    V_cycle = V_matrix(:,i);
    p_cycle = p_matrix(:,i);
    W_all(i) = trapz(V_cycle, p_cycle); % Numerical integration (trapezoidal) to find work per cycle
    m_fuel_cycle = m_fuel_matrix(:,i); % Mass fuel for cycle i
end
V_avg = sum(V_cycle)/Ncycles;
W_cumm = sum(W_all); %Cummulative function of the Work
m_fuel_cumm = sum(m_fuel_cycle);
m_per_cycle = m_fuel_cumm/Ncycles;
W_per_cycle = W_cumm/Ncycles; %Average work


%% Calculating derivatives
[dVdCa, smooth_dpdCa, smooth_p, V_avg, Ca, p, V, Ca_single, NCa] = calculatingDerivatives(Ca_matrix, p_matrix, V_matrix, Ncycles);

[smooth_P,dpdCa] = Pegging_dpdCa(smooth_p,NCa, Ca);
%rest of this was moved to "testing loop split"

%% aROHR selection data saving
if timing == chosen_aROHR
    Ca_saved = Ca;
    Gamma_at_angle_saved = Gamma_at_angle;
    smooth_p_saved = smooth_p;
    V_cycle_saved = V_cycle;
    dVdCa_saved = dVdCa;
    dpdCa_saved = dpdCa;
    answer = 1;
else
    answer = 0;
end

end


%% Mass flow
[W_per_cycle, m_per_cycle, V_avg, mass_flow_fuel_cycle, AFR] = calculateMassflow(Cyl, RPM, selectedFuel, FuelTable, V_cycle, m_fuel_cycle, W_all);

%% Plot PV diagram
[V_cycle] = plotPVDiagram(V_matrix, p_matrix, 1);

%% Plot all cycles
Pressure_Crankangle(Ca_matrix, p_matrix, bara, CaIVC, CaEVO, 10);


%% Plot average Pressure
plotAveragePressure(V_matrix, smooth_p, dm, bara, 10);

%% Calculates temperature in engine
[T, smooth_P, Ca_2to3, mfuel, mtot, Elcompfuel, Mi, LHV, Ymix, NSpS, T_curr, Yair, AF] = calculateTemperature(Ca, smooth_p, Cyl, SpS, FuelTable, V_cycle, CaIVC, Ca_single);

%% Gamma Calculation
[Gamma_at_angle] = calculateGamma(Ca_2to3, mfuel, mtot, Elcompfuel, SpS, Mi, FuelTable, LHV, Ymix, Runiv, T, NSpS, T_curr, Yair, Vmin, Vdiff, CaIVC, Cyl, Ca_single, V_cycle, smooth_P);

%% testing loop split

%The contents of this for loop was taken from under the "calculating
%derivatives" section

for T = 1:length(DataSetTimings)
    timing = DataSetTimings(T);
    
%Taking data out of emmision data table
EmissionLoadIndex = 50;
filteredRows = dataEmission(dataEmission.Load == EmissionLoadIndex & dataEmission.InjectionTiming == timing, :);
avgCO2 = mean(filteredRows.CO2);
avgNOx = mean(filteredRows.NOx);
VolumeEmission = calcEmissionVol(CaEVO, Cyl, smooth_P, Gamma_at_angle) ;        %Cylinder volume when exhaust valve opens
%Calculates the KPI for each loaded file and adds it to an array

[Efficiency_all, BSCO2_all, BSNOx_all, BSFC_all] = KPI_function(V_cycle, W_per_cycle,avgCO2,avgNOx, VolumeEmission,FuelTable,selectedFuel,smooth_p) ;
KPI_index_injection = timing;
Efficiency(T) = Efficiency_all;
BSCO2(T) = BSCO2_all;
BSNOx(T) = BSNOx_all;
BSFC(T) = BSFC_all;
injections(T) = DataSetTimings(T);


% save to struct or smthg data(i) 
end

%% Compute aROHR
if answer == 1 
    aROHR = computeAROHR(Ca_saved, Gamma_at_angle_saved, smooth_p_saved, V_cycle_saved, dVdCa_saved, dpdCa_saved);
else
    aROHR = computeAROHR(Ca, Gamma_at_angle, smooth_p, V_cycle, dVdCa, dpdCa);
end
%% Plot aROHR
f4 = figure(4);
set(f4, 'Position', [400 400 800 400]); % Figure size
plot(Ca(:, 1), aROHR, 'LineWidth', 1); % Plot averaged aROHR
xlabel('Crank Angle (°)');
ylabel('aROHR [J/°CA]');
legend()
xlim([-45 135]);
grid on;
title('Apparent Rate of Heat Release (aROHR) vs Crank Angle, changing gamma');

%% Find the Indices for CaSOI and CaEVO
[Ca_from_start, aROHR_from_start] = indices(Ca, CaSOI, CaEVO, aROHR);

%% Perform Cumulative Integration
aHR = cumtrapz(Ca_from_start, aROHR_from_start); % Cumulative heat release

%% Plot aHR
idx_BDC = NCa/4;
plotAHR(Ca_from_start, aHR, CaSOI, idx_BDC, p)

%% KPI Graphs

plotKPI(injections, Efficiency, BSCO2, BSNOx, BSFC);

%% Engine cycle
calculateEngineTemp(Cyl, SpS, Yair, AF, LHV, V_cycle)




