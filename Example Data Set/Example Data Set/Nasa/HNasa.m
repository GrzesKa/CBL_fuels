function [H,HU] = HNasa(T,Sp)
%HNasa(T,Sp):: computes internal specific enthalpy of species at Temp T
% 
%   Input: T, temperature
%          Sp, thermal database entry. So-called NASA polynomials
global Runiv
H=zeros(size(T));HU=H;
if (isempty(Runiv))
    fprintf('[HNasa] Assign global Runiv\n');
    return
end
if (isnan(Sp.Ts))
    Tl=T;
    a=Sp.Pol(1,:);
    HU=a(6) + a(1).*Tl+(a(2)/2).*Tl.^2+(a(3)/3).*Tl.^3+(a(4)/4).*Tl.^4+(a(5)/5).*Tl.^5;        % Formula 5.4 of lecture notes    
else
    ilow = (T <= Sp.Ts);
    Tl=T(ilow);
    a=Sp.Pol(1,:);
    HU(ilow)=a(6) + a(1).*Tl+(a(2)/2).*Tl.^2+(a(3)/3).*Tl.^3+(a(4)/4).*Tl.^4+(a(5)/5).*Tl.^5;        % Formula 5.4 of lecture notes
    
    ihigh = (T > Sp.Ts);
    Tl=T(ihigh);
    a=Sp.Pol(2,:);
    HU(ihigh)=a(6) + a(1).*Tl+(a(2)/2).*Tl.^2+(a(3)/3).*Tl.^3+(a(4)/4).*Tl.^4+(a(5)/5).*Tl.^5;        % Formula 5.4 of lecture notes
end
H=HU.*Runiv/Sp.Mass;
end

