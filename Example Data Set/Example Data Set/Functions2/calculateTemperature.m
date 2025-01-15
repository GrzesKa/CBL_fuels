function [T, Ca_2to3, mfuel, mtot, Elcompfuel, Mi, LHV, Ymix, NSpS, T_curr, Yair, AF, minimumsmooth_P] = calculateTemperature(Ca, smooth_P, Cyl, SpS, FuelTable, V_cycle, CaIVC, Ca_single)


min_value = 100000;

minimumsmooth_P = max (smooth_P, min_value);



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
MAir = Xair*Mi'*1000;                            %Molair mass of air 
Yair = Xair.*Mi/MAir;                       %Mass fraction of air
Xmix = Xair;
Ymix = Yair;
LHV = FuelTable.LHV(strcmp(FuelTable.Fuel, 'Diesel'));  %LHV of the fuel
Ca_2to3 = -3.2:0.2:149;                     %crank angle injection angle to exhaust valve opening
AF = 14.7;                                  %Air to fuel ratio
mfuel = mtot/length(Ca_2to3)/AF;            %Averaged fuel injection per increment
Elcompfuel = [SpS(1).Elcomp];               %Elemental composition of the fuel

% Corrected mass of air
mAir = 1.204 * Vbot * 1000; % Volume at BDC, mass in g

T_initial = 293; % Initial temperature [K]
R = 8.314; % Universal gas constant [J/(mol*K)]
T = zeros(size(Ca_single)); % Initialize temperature array
T(1) = T_initial; % Set initial temperature

% Corrected specific gas constant for air
R_air = R / MAir; % J/(kg*K)


CaIVC = -135;  %Intake valve closes
CaEVO = 149;   %Exhaust valve opens

% Corrected temperature calculation (only between CaIVC and CaEVO)
last_valid_temperature = NaN; % Initialize a variable to store the last valid temperature

for i = 1:length(Ca)
    
    if Ca(i) < CaIVC
        T(i) = 293; % Set temperature to 273 K before the intake valve closes
        continue;   % Skip the remaining calculations for this iteration
    end

    % Check if the crank angle is after the exhaust valve opens
    if Ca(i) > CaEVO
        T(i) = last_valid_temperature; % Set to the last calculated temperature
        continue;   % Skip the remaining calculations for this iteration
    end

    % Calculate current cylinder volume
    V_curr = CylinderVolume(Ca(i), Cyl); % Volume at current crank angle
    
    % Fetch current pressure
    P_curr = minimumsmooth_P(i); % Pressure at current crank angle

    % Check for valid volume
    if V_curr <= 0
        warning('Invalid volume at index %d: %f. Skipping...', i, V_curr);
        T(i) = NaN; % Set temperature to NaN for invalid volumes
        continue;
    end

    % Calculate temperature using the ideal gas law
    T_curr = P_curr * V_curr / (R * (mAir / 28.9647)); % [K]
    T(i) = T_curr; % Store temperature
    last_valid_temperature = T_curr; % Update the last valid temperature
end

% Plotting the temperature vs crank angle
figure(6);
plot(Ca, T, 'LineWidth', 1.5);
xlabel('Crank Angle (°)');
ylabel('Temperature [K]');
xlim([-180 180]);
title('Cylinder Temperature vs Crank Angle (Last Valid Temperature After Exhaust)');
grid on;
