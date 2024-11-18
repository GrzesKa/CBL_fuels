function [Cp,CpU] = CpNasa(T,Sp)
%CpNasa(T,Sp):: computes heat capacity of species at Temp T
% 
%   Input: T, temperature
%          Sp, thermal database entry. So-called NASA polynomials
global Runiv
Cp=zeros(size(T));CpU=Cp;
if (isempty(Runiv))
    fprintf('[CpNasa] Assign global Runiv\n');
    return
end
if (isnan(Sp.Ts))
    Tl=T;
    a=Sp.Pol(1,:);
    CpU=a(1)+a(2).*Tl+a(3).*Tl.^2+a(4).*Tl.^3+a(5).*Tl.^4;        % Formula 5.4 of lecture notes
else
    ilow = (T <= Sp.Ts);
    Tl=T(ilow);
    a=Sp.Pol(1,:);
    CpU(ilow)=a(1)+a(2).*Tl+a(3).*Tl.^2+a(4).*Tl.^3+a(5).*Tl.^4;        % Formula 5.4 of lecture notes
    
    ihigh = (T > Sp.Ts);
    Tl=T(ihigh);
    a=Sp.Pol(2,:);
    CpU(ihigh)=a(1)+a(2).*Tl+a(3).*Tl.^2+a(4).*Tl.^3+a(5).*Tl.^4;        % Formula 5.4 of lecture notes
end
Cp=CpU.*Runiv/Sp.Mass;
end

