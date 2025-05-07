function ds = newton3d_J2_moon(t, s)
%NEWTON3D_J2 Newtonian gravity trajectory propagation with J2 and lunar perturbations
%   Detailed explanation goes here

% values hardcoded for Earth:
mu = 3.9860043550702260E+14; % m^3/s^2
J2 = 1.08263E-3; 
Re = 6378100; % m

% values hardcoded for the Moon:
mu_moon = 4.9048695E12; % m^3/s^2
% hard coded based off of poor implementation:
jd_0 = 2457338.500000000;

jd_now = jd_0 + t / (3600 * 24);
mjd_now = jd_now - 2400000.5;

r = s(1:3);
r_mag = norm(r);
x = r(1);
y = r(2);
z = r(3);

v = s(4:6);
a_grav = -mu * r/(r_mag^3);


factor1 = 1.5 * J2 * (Re/r_mag)^2;
factor2 = 5*(z^2/r_mag^2);

% from page 79 of Alfriend
xdd_j2 = (mu * x / r_mag^3) * (factor1 * (factor2 - 1));
ydd_j2 = (mu * y / r_mag^3) * (factor1 * (factor2 - 1));
zdd_j2 = (mu * z / r_mag^3) * (factor1 * (factor2 - 3));

a_J2 = [xdd_j2; ydd_j2; zdd_j2];

r_moon = Moon(mjd_now);
r_sat_moon = r_moon - r;
r_sat_moon_mag = norm(r_sat_moon);
a_moon = mu_moon * r_sat_moon/(r_sat_moon_mag^3);

% a_drag = drag(r, v);

ddr = a_grav + a_J2 + a_moon;
ds = [v; ddr];
end