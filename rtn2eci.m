function [r_eci] = rtn2eci(r_chief_eci, v_chief_eci, rho_rtn)
%RTN2ECI Computes ECI position of a vector in RTN

r0 = r_chief_eci(:);
v0 = v_chief_eci(:);

rhat = r0/norm(r0);
h0 = cross(r0, v0);
nhat = h0/norm(h0);
that = cross(nhat, rhat);
C_RTN_to_ECI = [rhat, that, nhat];

rho_eci = C_RTN_to_ECI * rho_rtn;

r_eci = r_chief_eci + rho_eci;
end