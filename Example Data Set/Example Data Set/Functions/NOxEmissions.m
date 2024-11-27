function bsNOx = NOxEmissions(m_dot, P, selectedFuel, Fueltable)
%P =        load on engine
%m_dot =    mass flow rate of the fuel 

%Find NOx emmisions per kg of fuel value from table
NOx = FuelTable.NOx(strcmp(FuelTable.Fuel, selectedFuel));

m_dotNOx = m_dot*NOx; %mass flow of CO2 out of the engine
bsNOx = m_dotNOx/P;   %Brake specific emmision of CO2

end