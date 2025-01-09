function VolumeEmission = calcEmissionVol(CaEVO, Cyl, smooth_P)

%valve opens at crank angle 149
V_EVO = CylinderVolume(CaEVO, Cyl);
P_EVO = smooth_P(2545); 
Pamb = 101325;
gamma_EVO = 1.34 %Gamma_at_angle(2545,2);
%testing push
VolumeEmission = V_EVO * (P_EVO / Pamb)^(1/gamma_EVO);

end