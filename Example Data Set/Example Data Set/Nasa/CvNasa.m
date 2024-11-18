function [Cv,CvU] = CvNasa(T,Sp)
%CvNasa(T,Sp):: computes heat capacity of species at Temp T
% 
%   Input: T, temperature
%          Sp, thermal database entry. So-called NASA polynomials
global Runiv
if (isempty(Runiv))
    fprintf('[CvNasa] Assign global Runiv\n');
    return
end
[~,CpU]=CpNasa(T,Sp);
CvU = CpU-1;
Cv=CvU.*Runiv/Sp.Mass;
end

