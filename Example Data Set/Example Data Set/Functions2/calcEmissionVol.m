function VolumeEmission = calcEmissionVol(CaEVO, Cyl, smooth_P, Gamma_at_angle)

%valve opens at crank angle 149
V_EVO = CylinderVolume(CaEVO, Cyl);
P_EVO = smooth_P(2545); 
Pamb = 101325;
gamma_EVO = Gamma_at_angle(2545,2); %1.34;
%testing push
VolumeEmission = V_EVO * (P_EVO / Pamb)^(1/gamma_EVO);

end