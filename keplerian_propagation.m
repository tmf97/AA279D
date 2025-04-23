function [nus] = keplerian_propagation(nu, e, a, t)
%KEPLERIAN PROPAGATION 
%   Propagates true anomaly
%   Inputs:
%     nu - true anomaly [rad]
%     e - eccentricity [-]
%     a - semi-major axis [m]
%     t - time steps array [s]
%   Outputs:
%     nus - 

MU = 3.9860043550702260E+14; % m^3/s^2
M = ecc2mean(true2ecc(nu, e),e);

% propagate orbit
n = mean_motion(MU, a);
Ms = n.*t + M;
nus = zeros(size(t));
for i = 1:length(t)
    E = mean2ecc(Ms(i), e, 1e-12);
    nu = ecc2true(E, e);
    nus(i) = nu;
end