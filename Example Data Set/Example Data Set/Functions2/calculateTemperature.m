function [T, smooth_P, Ca_2to3, mfuel, mtot, Elcompfuel, Mi, LHV, Ymix, NSpS, T_curr, Yair, AF] = calculateTemperature(Ca, smooth_p, Cyl, SpS, FuelTable, V_cycle, CaIVC, Ca_single)
desired_value = 1e5;
% Pegging at 1 bar at BDC
for i = 1:length(smooth_p)
    % Calculate the difference between the original first value and the adjusted BDC value
    diff = smooth_p(180) - desired_value;
    
    % Adjust the current element by adding the calculated difference
    smooth_P(i) = smooth_p(i) - diff;
end
smooth_P=smooth_P(:);


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

% Corrected temperature calculation
for i = 1:length(Ca)
    V_curr = CylinderVolume(Ca(i), Cyl); % Volume at current crank angle
    P_curr = smooth_P(i); % Pressure at current crank angle
    
    % Check for valid V_curr and P_curr
    if V_curr <= 0
        warning('Invalid volume at index %d: %f. Skipping...', i, V_curr);
        continue;
    end
   

    T_curr = P_curr * V_curr / (R_air * mAir); % Temperature using ideal gas law [pV=nRT]
    T(i) = T_curr; % Store temperature
end



% Plotting Temperature vs Crank Angle
figure(6);
plot(Ca, T, 'LineWidth', 1.5);
xlabel('Crank Angle (°)');
ylabel('Temperature [K]');
xlim([-180 180]);
title('Cylinder Temperature vs Crank Angle');
grid on;