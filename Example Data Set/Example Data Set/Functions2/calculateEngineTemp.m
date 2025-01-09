function calculateEngineTemp(Cyl, SpS, Yair, AF, LHV, V_cycle)
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