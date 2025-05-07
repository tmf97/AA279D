function Phi = STM_kep(aoes, ti, tf, mu)
%STM_kep Summary of this function goes here
%   Detailed explanation goes here
a = aoes(1);
e = aoes(2);
i = aoes(3);
Omega = aoes(4);
w = aoes(5);
nu = aoes(6);

n = mean_motion(mu, a);
dt = tf - ti;
dM = n*dt;

Phi = eye(6);
Phi(2,1) = -1.5*dM;
end