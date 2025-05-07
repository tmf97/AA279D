function dt = orbit_timing(aoes, nuf, mu)
%ORBIT_TIMING Gives the time elapsed between two true anomalies
% Only handles a max delta-true-anomaly of 2pi

a = aoes(1);
e = aoes(2);
% i = aoes(3);
% raan = aoes(4);
% w = aoes(5);
nu0 = aoes(6);

M0 = true2mean(nu0, e);
Mf = true2mean(nuf, e);
dM = Mf - M0;
dM = wrapTo2Pi(dM);

n = mean_motion(mu, a);

dt = dM/n;
end