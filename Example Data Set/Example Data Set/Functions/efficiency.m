function efficiency = efficiency(m_dot, P, selectedFuel, FuelTable)

RPM = 1500;%cycles per minute
%P =        load on engine
%m_dot =    mass flow rate of the fuel 

%Find lower heating value from table
LHV = FuelTable.LHV(strcmp(FuelTable.Fuel, selectedFuel));

s = 60/RPM;%seconds per cycle
W = P*s;   %Work per cycle
Efuel = m_dot*s*LHV; %Energy content of the fuel
efficiency = W/Efuel;
end