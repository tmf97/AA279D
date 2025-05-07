function a_drag = drag(r, v)
%DRAG Computes Drag force on an MMS spacecraft in orbit around Earth
%   Only valid for altitudes above 1000km.

% Exponential atmosphere
% Data from Wertz (1978)
h = norm(r) - 6378000;
h0 = 1000000;
rho0 =  3.019e-15;
H =  268000;
rho = rho0 * exp(-(h-h0)/H);

% min cross-secitional area from https://mms.gsfc.nasa.gov/spacecraft.html 
A = 3.5 * 1.2; 

% page 50 of https://ntrs.nasa.gov/api/citations/19720018349/downloads/19720018349.pdf
Cd = 10.5;

m =  1250; % kg, approx mass https://mms.gsfc.nasa.gov/spacecraft.html

we = 7.292115e-5; %rad/s rotation rate of the Earth

v_rel = v - cross([0;0;we], r);

a_drag = -0.5 * (Cd*A/m) * rho * v_rel * norm(v_rel);
end