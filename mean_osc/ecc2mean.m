% ECC2MEAN computes mean anomaly given eccentric anomaly and eccentricity.
%
% Inputs:
%   E - eccentric anomaly [rad]
%   e - eccentricity [-]
%
% Outputs:
%   M - mean anomaly [rad]

function M = ecc2mean(E, e)

E = mod(E, 2*pi);

% Kepler's equation
M = E - e*sin(E);

end
