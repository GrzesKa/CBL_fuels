function [S,SU] = SNasa(T,Sp)
%SNasa(T,Sp):: computes internal specific entropy of species at Temp T
% 
%   Input: T, temperature
%          Sp, thermal database entry. So-called NASA polynomials
global Runiv
S=zeros(size(T));SU=S;
if (isempty(Runiv))
    fprintf('[SNasa] Assign global Runiv\n');
    return
end
if (isnan(Sp.Ts))
    Tl=T;
    a=Sp.Pol(1,:);
    SU=  a(7) + a(1).*log(Tl) +(a(2)/1).*Tl.^1+(a(3)/2).*Tl.^2+(a(4)/3).*Tl.^3+(a(5)/4).*Tl.^4;        % Formula 5.4 of lecture notes
else
    ilow = (T <= Sp.Ts);
    Tl=T(ilow);
    a=Sp.Pol(1,:);
%     SU(ilow)=a(7) + a(1).*log(Tl) + (a(2)/1) +(a(3)/2).*Tl.^1+(a(4)/3).*Tl.^2+(a(5)/4).*Tl.^3;        % Formula 5.4 of lecture notes
    SU(ilow)= a(7) + a(1).*log(Tl) +(a(2)/1).*Tl.^1+(a(3)/2).*Tl.^2+(a(4)/3).*Tl.^3+(a(5)/4).*Tl.^4;        % Formula 5.4 of lecture notes
    
    ihigh = (T > Sp.Ts);
    Tl=T(ihigh);
    a=Sp.Pol(2,:);
%     SU(ihigh)=a(7) + a(1).*log(Tl) + (a(2)/1) +(a(3)/2).*Tl.^1+(a(4)/3).*Tl.^2+(a(5)/4).*Tl.^3;        % Formula 5.4 of lecture notes
    SU(ihigh)= a(7) + a(1).*log(Tl) +(a(2)/1).*Tl.^1+(a(3)/2).*Tl.^2+(a(4)/3).*Tl.^3+(a(5)/4).*Tl.^4;        % Formula 5.4 of lecture notes
end
S=SU.*Runiv/Sp.Mass;
end
% ist = 5;
% cH=a(ist)/(ist-1);
% for ll=(ist-1):-1:2,
%     cH=cH*T+a(ll)/(ll-1);
% end;
% ll=7;
% cH=cH+a(1)*log(T)+a(ll);
