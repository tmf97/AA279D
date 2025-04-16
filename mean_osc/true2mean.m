% TRUE2MEAN computes mean anomaly given true anomaly and eccentricity.
%
%   Inputs:
%     f - true anomaly [rad]
%     e - eccentricity [-]
%
%   Outputs:
%     M - mean anomaly [rad]

function M = true2mean(f, e)

E = true2ecc(f,e);
M = ecc2mean(E,e);

end
