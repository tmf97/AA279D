function [r_eci, v_eci] = koe2pv(koes, mu)
%KOE2PV Convert Keplerian Orbital Elements to a position/velocity in ECI

a = koes(1);
e = koes(2);
i = koes(3);
Omega = koes(4);
w = koes(5);
nu = koes(6);

E = true2ecc(nu, e); % eccentric anomaly
n = mean_motion(mu, a); % mean motion

% calculate pos/vel in perifocal frame using Kepler's equations
xp = a*(cos(E)-e);
b = a*sqrt(1-e^2);
yp = b * sin(E);
zp = 0;

Edot = n / (1-e*cos(E));

vxp = -a * Edot*sin(E);
vyp = b*Edot*cos(E);
vzp = 0;

rp = [xp; yp; zp];
vp = [vxp; vyp; vzp];
pci = pfCeci(Omega, i, w);
r_eci = pci * rp;
v_eci = pci * vp;

end