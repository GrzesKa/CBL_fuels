function VolumeEmission = calcEmissionVol(CaEVO, Cyl, smooth_P, gammaEVO)

%valve opens at crank angle 149
V_EVO = CylinderVolume(CaEVO, Cyl);
P_EVO = smooth_P(2545); 
Pamb = 101325;
%testing push
VolumeEmission = V_EVO * (P_EVO / Pamb)^(1/gammaEVO);

end