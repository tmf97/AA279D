function dv_lb = dv_lower_bound(aDdroe, a, e, dt, MU)
%DV_LOWER_BOUND Computes in-plane reconfiguration delta-V lower bound

aDda = aDdroe(1);
aDdlambda = aDdroe(2);
aDde = aDdroe(3:4);
n = mean_motion(MU, a);
eta = sqrt(1-e^2);
dM = n*dt;

Dda_effect = abs(aDda) / (2*(1+e));
Ddlambda_effect = abs(aDdlambda) / (3*(1+e)*dM);
Dde_effect = norm(aDde) / sqrt(3*e^4 - 7*e^2 + 4);
dv_lb = n*eta * max([Dda_effect, Ddlambda_effect, Dde_effect]);
end