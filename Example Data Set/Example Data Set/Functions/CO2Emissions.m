function bsCO2 = CO2Emissions(m_dot, P, selectedFuel, Fueltable)
%P =        load on engine
%m_dot =    mass flow rate of the fuel 

%Find CO2 emmisions per kg of fuel value from table
CO2 = FuelTable.CO2(strcmp(FuelTable.Fuel, selectedFuel));

m_dotCO2 = m_dot*CO2; %mass flow of CO2 out of the engine
bsCO2 = m_dotCO2/P;   %Brake specific emmision of CO2

end