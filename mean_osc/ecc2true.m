% ECC2TRUE computes true anomaly given eccentric anomaly and eccentricity.
%
%   Inputs:
%     E - eccentric anomaly [rad]
%     e - eccentricity [-]
%
%   Outputs:
%     f - mean anomaly [rad]

function f = ecc2true(E, e)

E = mod(E, 2*pi);

f = mod(atan2(sin(E)*sqrt(1-e^2),cos(E)-e), 2*pi);

end
