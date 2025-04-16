% EQUIOE2KOE converts a set of equinoctial orbital elements to Keplerian 
% orbital elements.
% 
%   Inputs:
%     equioe - vector of equinoctial orbital elements:
%              a   - semi-major axis [m]
%              Psi - mean longitude [rad]
%              tq1 - [-]
%              tq2 - [-]
%              p1  - [-]
%              p2  - [-]
% 
%   Outputs:
%     koe - vector of Keplerian orbital elements:
%           a - semi-major axis [m]
%           e - eccentricity [-]
%           i - inclination [rad]
%           O - right ascension of the ascending node [rad]
%           w - argument of periapsis [rad]
%           M - mean anomaly [rad]

function koe = equioe2koe(equioe)

a = equioe(1);
Psi = equioe(2);
tq1 = equioe(3);
tq2 = equioe(4);
p1 = equioe(5);
p2 = equioe(6);

Omega = atan2(p2,p1);
i = atan2(2 * sqrt(p1^2 + p2^2), 1 - p1^2 - p2^2);

wtilde = atan2(real(tq2),real(tq1));
e = sqrt(tq1^2+tq2^2);

w = wtilde-Omega;
f = Psi-wtilde;

E = atan2(real(sin(f) * sqrt(1 - e^2)),real( cos(f) + e));
M = E - e*sin(E);

w = mod(w,2*pi);
M = mod(real(M),2*pi);
if (abs(M-2*pi) < eps)
    M = 0;
end

koe = [a; e; i; Omega; w; M];

end
