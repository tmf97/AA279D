function roe = pv2roe(pv_chief, pv_deputy)
%PV2ROE Converts PVs of chief and deputy spacecraft into ROEs
%   Nuff said
MU = 3.986004415e14; % Earth gravitational parameter [m^3/s^2]

koes_chief = pv2koe(pv_chief, MU);
koes_deputy = pv2koe(pv_deputy, MU);

roe = quasi_nonsingular_roe(koes_chief, koes_deputy);

end