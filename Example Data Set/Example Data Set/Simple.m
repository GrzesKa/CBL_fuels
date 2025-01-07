%% Info
% Toerental: 1500 RPM
% SOA van 4.2º voor TDC
% Resolutie van 0.2º CA
% Data voor 69 cycles (maximale van de Smetec, de OGO gensets kunnen in principe “onbeperkt” aan)
% 

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



%% Preallocate results structure
numLoads = 3;
EXP = 1; %Choose which experiment you want to load
results = struct ( 'Work', [], 'BSFC', [], 'Efficiency', [], 'aROHR', [], 'CumulativeHR', []) ;

%% Emission data

FullName = fullfile(sprintf('Data/EXP%d/EXP%d',EXP, EXP), sprintf('EmissionDataDiesel.txt')); % Adjust file naming as needed
dataEmission = readtable(FullName);
dataEmission.Properties.VariableNames = {'Load', 'InjectionTiming', 'CO', 'CO2', 'HC', 'O2', 'NOx'};
%% Loop through all the data
for p = 2:4
for T = 4:10

timing = T*2;
% Emmision dataset sorting


% Load dataset
FullName = fullfile(sprintf('Data/EXP%d/EXP%d/T%d',EXP, EXP, timing), sprintf('P%dT%d.txt', p,timing)); % Adjust file naming as needed
dataIn = table2array(readtable(FullName));

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
    m_fuel_cycle = m_fuel_matrix(:,i); % Mass fuel for cycle i
end
V_avg = sum(V_cycle)/Ncycles;
W_cumm = sum(W_all); %Cummulative function of the Work
m_fuel_cumm = sum(m_fuel_cycle);
m_per_cycle = m_fuel_cumm/Ncycles;
W_per_cycle = W_cumm/Ncycles; %Average work

%Taking data out of emmision data table
EmissionLoadIndex = (p-1)*20+10;
filteredRows = dataEmission(dataEmission.Load == EmissionLoadIndex & dataEmission.InjectionTiming == timing, :);
avgCO2 = mean(filteredRows.CO2)
avgNOx = mean(filteredRows.NOx)
VolumeEmission = CylinderVolume(CaEVO,Cyl);         %Cylinder volume when exhaust valve opens
%Calculates the KPI for each loaded file and adds it to an array
[Efficiency_all, BSCO2_all, BSNOx_all, BSFC_all] = KPI_function(V_cycle, W_per_cycle,avgCO2,avgNOx, VolumeEmission);
KPI_index_injection = T-3;
KPI_index_load = p-1;
Efficiency(KPI_index_load,KPI_index_injection) = Efficiency_all;
BSCO2(KPI_index_load,KPI_index_injection) = BSCO2_all;
BSNOx(KPI_index_load,KPI_index_injection) = BSNOx_all;
BSFC(KPI_index_load,KPI_index_injection) = BSFC_all;
injections(KPI_index_injection) = T*2;
end
end

%Mass flow
intake_pressure = 100000; %Pa ambient pressure 1 bar
intake_temp = 293.15; %Kelvin ambient temperature
spec_gasct_air = 287; %J/mol*K
volume_displaced_cycle = (pi/4) * Cyl.Bore\mm * Cyl.Stroke\mm;
volumetric_flow = volume_displaced_cycle * 0.5 * (RPM/60);
density_air = intake_pressure/(intake_temp*spec_gasct_air);
mass_flow_air = volumetric_flow * density_air;

m_air = (mass_flow_air*2)/(RPM/60);
AFR = m_air./m_fuel_cycle;

mass_flow_fuel = mass_flow_air./AFR; %air to fuel ratio
mass_flow_fuel_cycle = mass_flow_fuel/12.5;


% Display or save results
disp('Work per cycle:');
disp(W_per_cycle);
disp('mpercycle:');
disp(m_per_cycle);
disp('Volume avg');
disp(V_avg);



% Plot p-V diagram for one cycle (e.g., cycle 1)
figure;
plot(V_matrix(:,1), p_matrix(:,1));
xlabel('Volume [m^3]');
ylabel('Pressure [Pa]');
title('p-V Diagram for Cycle 1');
grid on;

LHV = FuelTable.LHV(strcmp(FuelTable.Fuel, selectedFuel));
disp('LHV ::');
disp(LHV);

%[Efficiency_all, BSCO2_all, BSNOx_all, BSFC_all] = KPI_function(V_cycle, W_per_cycle);

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


%(function) for loop
%In: mfuel,mtot Ca, Xmix, Ymix,CaIVC, T, 
%Out: Xmix, Ymix, gamma, temperature, QLHV





%% Compute Volume and Derivatives

% Define variables
Ca = Ca_matrix;       % Use all crank angle data
p = p_matrix;         % Use all pressure data
V = V_matrix;         % Use all volume data
Ca_single = Ca(:, 1); % Use the crank angle array (same for all cycles)
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
    filtered_averaged_p(i) = filtered_averaged_p(i)-diff;
end
smooth_P=smooth_P(:);

Ca_single = Ca(:, 1);
dpdCa = zeros(size(Ca_single));

dpdCa = gradient(smooth_P, Ca_single);

dpdCa = sgolayfilt(dpdCa, 5, 37);

f6 = figure(6);
subplot(2,2,1)
hold on;
plot(Ca_single,dpdCa) %filtered_averaged_p/bara)
plot(Ca_single,filtered_averaged_dpdCa)
plot(Ca_single,smooth_dpdCa)
legend('dpdCa','filtered_averaged_dpdCa','smooth_dpdCa')
xlim([-45 135])
hold off;
subplot(2,2,2)
plot(Ca_single,V_cycle)
xlim([-45 135])
legend()
subplot(2,2,3)
hold on;
plot(Ca_single,filtered_averaged_p)
plot(Ca_single, smooth_P)
plot(Ca_single, smooth_p)
legend('filtered_averaged_p','smooth_P','smooth_p')
xlim([-45 135])
hold off;
subplot(2,2,4)
hold on;
plot(Ca_single,dpdCa)

xlim([-45 135])
legend()
hold off;

Vmin = min(V_cycle(:)); %Find the smallest value from the Volume matrix
B   = Cyl.Bore; %Defined in the CylinderVolume function
S   = Cyl.Stroke; %Defined in the CylinderVolume function
Vdiff = pi / 4 * B^2 * S; % Displacement volume
NSpS = length(SpS); %calculates the length of the SpS matrix
Pamb = 101325; % pressure in Pa
B   = Cyl.Bore;
S   = Cyl.Stroke;


% known volumes 
Vbot = Vmin + Vdiff; % [m^3] the volume of the cylinder at BDC
Vtop = Vmin; %[m^3] the volume of the cylinder at TDC

VIntakeClose = CylinderVolume(CaIVC,Cyl);   %Volume of cylinder after intake valve closes
DensityAir = 1.204;                         %kg/m3 density of air at atmospheric pressure
MassAirIntake = VIntakeClose * DensityAir;  %Mass of the air after intake valve closes
mtot = MassAirIntake;
Xair = [0 0.21 0 0 0.79];                   %Mol fraction of air
Mi = [SpS.Mass];                            %Molair mass call from Nasa table
MAir = Xair*Mi';                            %Molair mass of air 
Yair = Xair.*Mi/MAir;                       %Mass fraction of air
Xmix = Xair;
Ymix = Yair;
LHV = FuelTable.LHV(strcmp(FuelTable.Fuel, 'Diesel'));  %LHV of the fuel
Ca_2to3 = -3.2:0.2:149;                     %crank angle injection angle to exhaust valve opening
AF = 14.7;                                  %Air to fuel ratio
mfuel = mtot/length(Ca_2to3)/AF;            %Averaged fuel injection per increment
Elcompfuel = [SpS(1).Elcomp];               %Elemental composition of the fuel

% Corrected mass of air
mAir = 1.293 * Vbot; % Volume at BDC, mass in g

T_initial = 293; % Initial temperature [K]
R = 8314; % Universal gas constant [J/(mol*K)]
T = zeros(size(Ca_single)); % Initialize temperature array
T(1) = T_initial; % Set initial temperature

% Corrected specific gas constant for air
R_air = R / MAir; % J/(kg*K)

% Corrected temperature calculation
Gamma_at_angle = zeros(length(Ca_single), 2); %Initiates matrix of length(Ca_2to3 by 2

for i = 1:length(Ca_single)
    Gamma_at_angle(i,1) = Ca_single(i);           %Appends each crank angle to the matrix
    V_curr = V_cycle(i); % Volume at current crank angle
    P_curr = smooth_P(i); % Pressure at current crank angle
    
    % Check for valid V_curr and P_curr
    if V_curr <= 0
        warning('Invalid volume at index %d: %f. Skipping...', i, V_curr);
        continue;
    end
    


    

    T_curr = P_curr * V_curr / (R_air * mAir); % Temperature using ideal gas law [pV=nRT]
    T(i) = T_curr; % Store temperature
    % Calculate Cp and Cv for the current temperature (Sander's new code)
    for j = 1:NSpS
        Cp_T(j) = CpNasa(T(i), SpS(j));
        Cv_T(j) = CvNasa(T(i), SpS(j)); 
    end
    
    % Weighted averages based on air composition
    Cp_avg = Yair * Cp_T';
    Cv_avg = Yair * Cv_T';
    R_air = Cp_avg-Cv_avg;
    gamma_curr = Cp_avg/Cv_avg;
    Gamma_at_angle(i,2) = gamma_curr;    %Appends the gamma value in the second column
end





% Plotting Temperature vs Crank Angle
figure;
plot(Ca, T, 'LineWidth', 1.5);
xlabel('Crank Angle (°)');
ylabel('Temperature [K]');
xlim([-360 360]);
title('Cylinder Temperature vs Crank Angle');
grid on;


% Plotting Gamma vs Crank Angle
figure;
plot(Ca_single, Gamma_at_angle(:,2), 'LineWidth', 1.5);
xlabel('Crank Angle (°)');
ylabel('Gamma [-]');
xlim([-360 360]);
title('Gamma value vs Crank Angle');
grid on;

%% Gamma calculations
%T = 800;

%For loop that incrementally adds the fuel every .2 crank angle increments
%and calculates the new composition, gamma value and temperature
%Gamma_at_angle = zeros(length(Ca_2to3), 2); %Initiates matrix of length(Ca_2to3 by 2

for i=1:length(Ca_2to3)
%Gamma_at_angle(i,1) = Ca_2to3(i);           %Appends each crank angle to the matrix
end

%Current mol values in the cylinder
nFueladd = mfuel / Mi(1);             % Mols of fuel
nO2current = mtot * Ymix(2) / Mi(2);  % Mols of O2               
nN2current = mtot * Ymix(3) / Mi(3);  % Mols of N2 (remains unchanged)
nCO2current = mtot * Ymix(4) / Mi(4); % Mols of CO2
nH2Ocurrent = mtot * Ymix(5) / Mi(5); % Mols of H2O
nTotalCurrent = nFueladd + nO2current + nN2current + nCO2current + nH2Ocurrent;
mtot = mtot + mfuel;                  %Total mass in the cylinder

%Mol going out                                    
nO2new = nO2current - (Elcompfuel(3)+Elcompfuel(2)/4)*nFueladd;    % new Mols of O2   
nN2new = nN2current;                                               % new Mols of N2 (remains unchanged)
nCO2new = nCO2current + nFueladd * Elcompfuel(3);                  % new Mols of CO2
nH2Onew = nH2Ocurrent + nFueladd * Elcompfuel(2)/2;                % new Mols of H2O
nNew = nO2new + nN2new + nCO2new + nH2Onew;                        


%Molar mass fractions of the products
Xmix(1) = 0;
Xmix(2) = nO2new/nNew;
Xmix(3) = nN2new/nNew;
Xmix(4) = nCO2new/nNew;
Xmix(5) = nH2Onew/nNew;

Mmix = Xmix * Mi';         % Mean molar mass of the mixture
Ymix = Xmix .* Mi / Mmix;  % Mass fractions of products
Rgmix = Runiv / Mmix;      % Specific gas constant of the mixture

% Alexander's code
%For loop calculates the Cp and Cv values for each element in the Array
%     for j=1:5                                                                                       
%         Cpi(:,j) = CpNasa(T,SpS(j));
%         Cvi(:,j) = CvNasa(T,SpS(j)); 
%     end
% Cp = Ymix*Cpi';                 %Calculates the Cp of the fuel using each component Cpi
% Cv = Ymix*Cvi';                 %Calculates the Cv of the fuel using each component Cvi
% Cpcheck = Cv + Rgmix;           %Checks if the  Cp value is consistent with the formula Cp = Cv + Rspecif
% Gamma = Cp/Cv;                  %Calculates gamma using Cp and Cv
% Gamma_at_angle(i,2) = Gamma;    %Appends the gamma value in the second column
% deltaT = LHV*mfuel/(Cp*mtot);   %Calculates
% T = T + deltaT;
% end

% Calculate Cp and Cv for the current temperature (Sander's new code)
    for j = 1:NSpS
        Cp_T(j) = CpNasa(T_curr, SpS(j));
        Cv_T(j) = CvNasa(T_curr, SpS(j)); 
    end
    
    % Weighted averages based on air composition
    Cp_avg = Yair * Cp_T';
    Cv_avg = Yair * Cv_T';
    gamma_curr = Cp_avg/Cv_avg;

%emissions
CO2_endmass = Ymix(4)*mtot;
disp('CO2 endmass in grams');
disp(CO2_endmass*1000);


%% Compute aROHR

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

% % Initialize aROHR

aROHR = zeros(1, length(Ca_single)); % Preallocate for speed, assuming numDatasets is defined
% Compute aROHR using the interpolated gamma
aROHR = (gamma_full ./ (gamma_full - 1)) .* smooth_p .* dVdCa(:,1) + ...
        (1 ./ (gamma_full - 1)) .* V_cycle .*dpdCa;
%% Plot aROHR

disp('aROHR size is');
disp(size(aROHR));

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
f3 = figure(3);
figure;
hold on; % Hold on to plot all datasets on the same figure
for datasetIndex = 1:numLoads
    plot(Ca_single, aROHR, 'DisplayName', 'Dataset 1');
    plot(Ca_single, aROHR, 'DisplayName', sprintf('Dataset %d', datasetIndex));
    hold on; % Hold on to plot all datasets on the same figure
end

xlabel('Crank Angle (°)');
ylabel('aROHR [J/°CA]');
xlim([-45 135]);
grid on;
title('Apparent Rate of Heat Release (aROHR) vs Crank Angle');
legend show; % Show legend to differentiate datasets
hold off;



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
xlim([CaSOI - 10, 180]); % Adjust x-axis limits
ylim([min(aHR) - 100, max(aHR) + 100]); % Adjust y-axis limits
grid on;
hold off;



%% KPI Graphs
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
title('Efficiency at Various Loads');

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


title('BSCO2 at Various Loads');

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
title('BSNox at Various Loads');

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
title('BSFC at Various Loads');

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

function [Efficiency_all, BSCO2_all, BSNOx_all, BSFC_all] = KPI_function(V_cycle, W_per_cycle,CO2, NOx,VolumeEmission)

RPM = 1500 ; % Rotations Per Minute of engine
Power_engine = ((W_per_cycle/1000) * RPM)/120 ; % Power of engine, work hardcoded due to errors
LHV_B7 = 43e6; % LHV of B7 diesel in J/kg
V_air =  max(V_cycle); % volume of air in 0.5 cycle * nr of cycles in 1 second
mass_fuel = 0.0001176; % mass of fuel being injected (data from sensors)
ro_air = 1.293 ; % Density air 
V_air_ps = V_air*RPM/120 ; % Volume of air per second
m_air = V_air_ps * ro_air ; % mass of air entering, calculated using v of air in 1 second * density
exhaust_massflow = mass_fuel + m_air; % total mass flowing
exhaust_massflow_grams = exhaust_massflow * 1000; % used for the brake specific KPI's in g/
Efuel = mass_fuel*LHV_B7 ; 

%% KPIs calculation 

Efficiency_all = Power_engine*1000/Efuel ; % work per cycle/Efuel (work in J)

%Molar mass
NOx_NO2 = 0.1; % Percentage of NOx that transforms into NO2
molar_CO2 = 44;   % g/mol
molar_NOx = (30 * (1 - NOx_NO2) + 46 * NOx_NO2); % Weighted molar mass of NOx (g/mol)
molar_H2O = 18;   % g/mol
molar_Ar = 40;
molar_N2 = 28;    % g/mol
molar_O2 = 32;    % g/mol


% Combustion volume fractions of diesel
emission_N2 = 0.76;         % volume fraction for N2
emission_O2 = 0.135;        % volume fraction for O2
emission_CO2 = 0.0525;      % volume fraction for CO2
emission_H20 = 0.05;        % volume fraction for H2O
emission_Ar = 0.008;
emission_NOx = 0.00125;     % volume fraction for NOx


% Total molar mass contribution
total_molar_mass = (emission_N2 * molar_N2) + (emission_O2 * molar_O2) + (emission_CO2 * molar_CO2) + (emission_H20 * molar_H2O) + (emission_Ar * molar_Ar) + (emission_NOx * molar_NOx);


% Volume-to-mass fraction conversion
mass_fraction_CO2 = (emission_CO2 * molar_CO2) / total_molar_mass;
mass_fraction_NOx = (emission_NOx * molar_NOx) / total_molar_mass;

% Mass flow rates OLD
%CO2_massflow = exhaust_massflow_grams * mass_fraction_CO2; 
%NOx_massflow = exhaust_massflow_grams * mass_fraction_NOx; 

% Mass flow rates NEW
DensityCO2 = 1.98e3;                                %Density of CO2 at 20 degrees celcius (g/m^3)

CO2_massflow = CO2/100*VolumeEmission*DensityCO2/100;
NOx_massflow = NOx*46.01/24*VolumeEmission/100;

%KPIs - x3600 to get to kWhr
BSCO2_all = (CO2_massflow*3600)/Power_engine;
BSNOx_all = (NOx_massflow*3600)/Power_engine;
BSFC_all = (mass_fuel*3600)/Power_engine;

disp(['Efficiency_all:', num2str(Efficiency_all)]);

disp(['BSCO2_all:', num2str(BSCO2_all), 'g/kWhr']);

disp(['BSNOx_all:', num2str(BSNOx_all), 'g/kWhr']);

disp(['BSFC_all:', num2str(BSFC_all), 'g/kWhr']);


end







%% Engine cycle

%variable definitions 
Vmin = min(V_cycle(:)); %Find the smallest value from the Volume matrix
B   = Cyl.Bore; %Defined in the CylinderVolume function
S   = Cyl.Stroke; %Defined in the CylinderVolume function
Vdiff = pi / 4 * B^2 * S; % Displacement volume
NSpS = length(SpS); %calculates the length of the SpS matrix
Pamb = 101325; % pressure in Pa

%The function itself
Vmin = min(V_cycle(:));
B   = Cyl.Bore;
S   = Cyl.Stroke;
Vdiff = pi / 4 * B^2 * S; % Displacement volume
NSpS = length(SpS);
Pamb = 101325; % pressure in Pa
V1 = Vmin + Vdiff; % [m^3] the volume of the cylinder at BDC
V2 = Vmin; %[m^3] the volume of the cylinder at TDC
T1 = 293; %  assumption room temperature
for i=1:NSpS
    Cp_T1(i) = CpNasa(T1,SpS(i));
    Cv_T1(i) = CvNasa(T1,SpS(i)); 
end

Cp_1 = Yair*Cp_T1'; % [J/(kg*K)]specific heat for constant pressure at temperature 1
Cv_1 = Yair*Cv_T1'; % [J/(kg*K)] specific heat for constant volume at temperature 1
gamma_T1 = Cp_1/Cv_1; % gamma from specific heat at temperature 1
R = Cp_1-Cv_1; % specific gas constant

% compression from step 1 to 2
mass_air = Pamb*V1/(T1*R); % [kg]the mass in the cylinder
Ca_1_2 = linspace(180, 360, 360); %crank angle 180 to 360 degrees (compression stroke)

% Calculate volumes at each crank angle
V1_2 = arrayfun(@(ca) CylinderVolume(ca, Cyl), Ca_1_2);
P1_2 = Pamb * (V1 ./ V1_2).^gamma_T1; % pressures using Poisson
P2 = P1_2(end); % [bar] final value of the compression stroke
T2 = P2*Vmin/(R*mass_air); % [K] temperature at point 2

% combustion from step 2 to 3
P3 = P2; %[bar] isobaric condition
m_fuel = mass_air/AF; % [kg] mass of fuel
Yfuel_mix = 1/(AF+1); % mass fraction fuel mixture
Yair_mix = AF/(AF+1); % mass fraction air mixture
Ymix = Yfuel_mix*[1 0 0 0 0] + Yair_mix*Yair;
for i=1:NSpS
    Cp_T2(i) = CpNasa(T2,SpS(i));
    
end
Cp2 = Ymix*Cp_T2'; % [J/(kg*K)]specific heat for constant pressure at temperature 1
Qin = m_fuel*LHV; % [J] heat into the system step 2 to 3
T3 = Qin/(Cp2.*(m_fuel+mass_air)) + T2; % 
V3 = V2*T3/T2; % using ideal gas law
for i=1:NSpS
    Cp_T3(i) = CpNasa(T3,SpS(i));
    Cv_T3(i) = CvNasa(T3,SpS(i)); 
end
Cp3 = Yair*Cp_T3'; % [J/(kg*K)]specific heat for constant pressure at temperature 3
Cv3 = Yair*Cv_T3'; % [J/(kg*K)] specific heat for constant volume at temperature 3
gamma_T3 = Cp3/Cv3; % gamma from specific heat at temperature 3
T4 = T3*(V2/V1)^(gamma_T3-1);

T1 %outputs temperatures in command window
T2
T3
T4

