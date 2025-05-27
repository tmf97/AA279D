function [rho_RTN, rho_dot_RTN_in_RTN] = pvs2rtn(r0, v0, r1, v1)
%PVS2RTN Convert Position/Velocity of Chief and Deputy to RTN
%   Compute position and velocity of r1/v1 in the reference frame of 0.
r0 = r0(:);
v0 = v0(:);
r1 = r1(:);
v1 = v1(:);

rho_ECI = r1-r0;
% define RTN frame
rhat = r0/norm(r0);
h0 = cross(r0, v0);
nhat = h0/norm(h0);
that = cross(nhat, rhat);
C_ECI_to_RTN = [rhat, that, nhat].';

% angular speed of the RTN frame in RTN
% fdot is the time derivative of true anomaly.
% r^2 * fdot = h
fdot = norm(h0) / (norm(r0)^2);

omega_RTN_in_ECI_in_RTN = [0;0; fdot];

rho_RTN = C_ECI_to_RTN * rho_ECI;

rho_dot_ECI = v1 - v0;

omega_ECI_in_RTN_in_ECI = C_ECI_to_RTN.' * -omega_RTN_in_ECI_in_RTN;

rho_dot_RTN_in_ECI = rho_dot_ECI + cross(omega_ECI_in_RTN_in_ECI, rho_ECI);
rho_dot_RTN_in_RTN = C_ECI_to_RTN * rho_dot_RTN_in_ECI;
end