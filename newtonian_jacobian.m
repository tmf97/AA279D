function [A] = newtonian_jacobian(r_vec, mu)
%NEWTONIAN_JACOBIAN Computes Jacobian of FODE
%   
% Inputs:
%   r_vec - 3x1 position vector [x; y; z] in m
%   mu    - Gravitational parameter (m^3/s^2)
%
% Output:
%   A     - 6x6 Jacobian matrix

r = norm(r_vec);
I3 = eye(3);
rrT = r_vec * r_vec.';  % Outer product r * r^T

dadr = mu * ( (-1 / r^3) * I3 + (3 / r^5) * rrT );

A = zeros(6,6);
A(1:3, 4:6) = I3;
A(4:6, 1:3) = dadr;
end