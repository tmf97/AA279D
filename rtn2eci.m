function [vec_eci] = rtn2eci(r_chief_eci, v_chief_eci, vec_rtn)
%RTN2ECI Gets ECI representation of rho_rtn vector

r0 = r_chief_eci(:);
v0 = v_chief_eci(:);

rhat = r0/norm(r0);
h0 = cross(r0, v0);
nhat = h0/norm(h0);
that = cross(nhat, rhat);
C_RTN_to_ECI = [rhat, that, nhat];

vec_eci = C_RTN_to_ECI * vec_rtn;

end