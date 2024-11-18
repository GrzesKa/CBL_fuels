function [U,UU] = UNasa(T,Sp)
%CvNasa(T,Sp):: computes heat capacity of species at Temp T
% 
%   Input: T, temperature
%          Sp, thermal database entry. So-called NASA polynomials
global Runiv
if (isempty(Runiv))
    fprintf('[UNasa] Assign global Runiv\n');
    return
end
[~,HU]=HNasa(T,Sp);
UU = HU-T;
U=UU.*Runiv/Sp.Mass;
end

