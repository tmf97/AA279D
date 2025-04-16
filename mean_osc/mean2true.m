% MEAN2TRUE computes true anomaly from eccentricity and mean anomaly. Note
% that this function relies on a Newton-Raphson method to numerically
% compute a value for true anomaly.
%
%   Inputs:
%     M   - mean anomaly [rad]
%     e   - eccentricity of orbit [-]
%     tol - tolerance for Newton-Raphson iterator
%
%   Outputs:
%     f - true anomaly [rad]

function f = mean2true(M, e, tol)

if nargin == 2
    tol = 10e-10;
end

E = mean2ecc(M,e,tol);
f = ecc2true(E,e);

end
