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
%% Load NASA maybe you need it at some point?
% Global (for the Nasa database in case you wish to use it).
global Runiv
Runiv = 8.314;
[SpS,El]        = myload('Nasa\NasaThermalDatabase.mat',{'Diesel','O2','N2','CO2','H2O'});
%% Engine geom data (check if these are correct)
Cyl.Bore                = 104*mm;
Cyl.Stroke              = 85*mm;
Cyl.CompressionRatio    = 21.5;
Cyl.ConRod              = 136.5*mm;
Cyl.TDCangle            = 180;
% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).
CaIVO = -355; %Intake valve opens
CaIVC = -135; %Intake valve closes
CaEVO = 149; %Exhaust valve opens
CaEVC = -344; %Exhaust valve closes
CaSOI = -3.2; %Start of Injection - CHANGE IF WE PLAY AROUND WITH IT IN THE EXPERIMENT
% Write a function [V] = CylinderVolume(Ca,Cyl) that will give you Volume
% for the given Cyl geometry. If you can do that you can create pV-diagrams
%% Load data (if txt file)
FullName        = fullfile('Data','ExampleDataSet.txt');
dataIn          = table2array(readtable(FullName));
[Nrows,Ncols]   = size(dataIn);                    % Determine size of array
NdatapointsperCycle = 720/0.2;                     % Nrows is a multitude of NdatapointsperCycle
Ncycles         = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
Ca              = reshape(dataIn(:,1),[],Ncycles); % Both p and Ca are now matrices of size (NCa,Ncycles)
p               = reshape(dataIn(:,2),[],Ncycles)*bara; % type 'help reshape' in the command window if you want to know what it does (reshape is a Matlab buit-in command
%% Plotting 
f1=figure(1);
set(f1,'Position',[ 200 800 1200 400]);             % Just a size I like. Your choice
pp = plot(Ca,p/bara,'LineWidth',1);                 % Plots the whole matrix
xlabel('Ca');ylabel('p [bar]');                     % Always add axis labels
xlim([-360 360]);ylim([0 50]);                      % Matter of taste
iselect = 10;                                       % Plot cycle 10 again in the same plot to emphasize it. Just to show how to access individual cycles.
line(Ca(:,iselect),p(:,iselect)/bara,'LineWidth',2,'Color','r');
YLIM = ylim;
% Add some extras to the plot
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b'); % Plot a vertical line at IVC. Just for reference not a particular reason.
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r'); % Plot a vertical line at EVO. Just for reference not a particular reason.
set(gca,'XTick',[-360:60:360],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('All cycles in one plot.')
%% pV-diagram
V = CylinderVolume(Ca(:,iselect),Cyl);
f2 = figure(2);
set(f2,'Position',[ 200 400 600 800]);              % Just a size I like. Your choice
subplot(2,1,1)
plot(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0 0.8]);ylim([0.5 50]);                      % Matter of taste
set(gca,'XTick',[0:0.1:0.8],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram'})
subplot(2,1,2)
loglog(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0.02 0.8]);ylim([0 50]);                      % Matter of taste
set(gca,'XTick',[0.02 0.05 0.1 0.2 0.5 0.8],...
    'YTick',[0.5 1 2 5 10 20 50],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram'})

%% Constants
%Cp = CpNasa(T, )
%Cv = CvNasa(T, )

gamma = 1.3;        % Ratio of specific heats, this is used for now

%% Compute Volume and Derivatives
V = CylinderVolume(Ca, Cyl); % Calculate volume at each crank angle

dVdCA = smoothdata(gradient(V, diff(Ca(1:2))), 'movmean', 5); % Approximation of dV/dCA

dpdCA = smoothdata(gradient(p, diff(Ca(1:2))), 'movmean', 5); % Approximation of dp/dCA


%% Compute aROHR
aROHR = -(gamma / (gamma - 1)) * p .* dVdCA - (1 / (gamma - 1)) * V .* dpdCA;


%% Plot aROHR
f3 = figure(3);
set(f3, 'Position', [400 400 800 400]); % Figure size
plot(Ca(:, iselect), aROHR(:, iselect), 'LineWidth', 1); % Plot for selected cycle
xlabel('Crank Angle (°)');
ylabel('aROHR [J/°CA]');
xlim([-45 135]);
%ylim([-20 100])
grid on;
title('Apparent Rate of Heat Release (aROHR) vs Crank Angle');

%% Find the Index for CaSOI
[~, idx_start] = min(abs(Ca(:, iselect) - CaSOI)); % Find index closest to CaSOI

% Slice data starting at CaSOI to the last data point
Ca_from_start = Ca(idx_start:end, iselect); % Crank angle from CaSOI to the end
aROHR_from_start = aROHR(idx_start:end, iselect); % aROHR from CaSOI to the end

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
xlim([-10 300]); % Adjust x-axis limits
ylim([-100 5000]); % Adjust y-axis limits
grid on;
hold off;
