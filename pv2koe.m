function koe = pv2koe(pv, mu)
%PV2KOE Convert position/velocity to Keplerian Orbital Elements
%   Detailed explanation goes here
r_eci = pv(1:3);
v_eci = pv(4:6);
r = norm(r_eci);
rhat = r_eci / r;
v = norm(v_eci);

h = cross(r_eci, v_eci);
hnorm = norm(h);

e_vec = cross(v_eci, h) / mu - rhat;
e = norm(e_vec);

n_vec = cross([0, 0, 1] , h);
n = norm(n_vec);

specific_energy = v^2/2 - mu/r;

a = -mu / (2 * specific_energy);

i = acos(h(3)/hnorm);
Omega = acos(n_vec(1)/n);
w = acos(dot(n_vec, e_vec)/(n*e));
nu = acos(dot(e_vec, r_eci)/(e*r));


if n_vec(2) < 0
    Omega = 2*pi - Omega;
end

if e_vec(3) < 0
    w = 2*pi - w;
end

if dot(r_eci, v_eci) < 0
    nu = 2*pi - nu;
end

koe = [a, e, i, Omega, w, nu];
end