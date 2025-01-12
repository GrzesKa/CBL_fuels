function [Gamma_at_angle] = calculateGamma(Ca_2to3, mfuel, mtot, Elcompfuel, SpS, Mi, FuelTable, LHV, Ymix, Runiv, T, NSpS, T_curr, Yair, Vmin, Vdiff, CaIVC, Cyl, Ca_single, V_cycle, smooth_P)
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
    


    
    %T_initial = 293 
    %T_curr = T_initial + (P_curr * V_curr / (R_air * mAir)); % Temperature using ideal gas law [pV=nRT]
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
    Gamma_at_angle(i,2) = gamma_curr;    %Appends the gamma value in the second column
end