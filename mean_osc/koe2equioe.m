% KOE2EQUIOE converts a set of Keplerian orbital elements to equinoctial orbital 
% elements.
% 
%   Inputs:
%     koe - vector of Keplerian orbital elements:
%           a - semi-major axis [m]
%           e - eccentricity [-]
%           i - inclination [rad]
%           O - right ascension of the ascending node [rad]
%           w - argument of periapsis [rad]
%           M - mean anomaly [rad]
% 
%   Outputs:
%     equioe - vector of equinoctial orbital elements:
%              a   - semi-major axis [m]
%              Psi - mean longitude [rad]
%              tq1 - [-]
%              tq2 - [-]
%              p1  - [-]
%              p2  - [-]

function equioe = koe2equioe(koe)

a = koe(1);
e = koe(2);
i = koe(3);
Omega = koe(4);
w = koe(5);
M = koe(6);

f = mean2true(M,e);
w_tilde = Omega + w;
Psi = w_tilde + f;
tq1 = e * cos(w_tilde);
tq2 = e * sin(w_tilde);
p1 = tan(i / 2) * cos(Omega);
p2 = tan(i / 2) * sin(Omega);

Psi = mod(Psi,2*pi);
if Psi > pi
   Psi = Psi-2*pi; 
end

equioe = [ a; Psi; tq1; tq2; p1; p2 ];

end
