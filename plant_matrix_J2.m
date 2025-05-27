function A = plant_matrix_J2(aoes, mu)
% PLANT_MATRIX_J2 
% computes the J2 + Keplerian plant matrix from eqn 23 of Koenig 2016

J2 = 1.08263E-3; 
Re = 6378100; % m


a = aoes(1);
e = aoes(2);
i = aoes(3);
% Omega = aoes(4);
w = aoes(5);
% nu = aoes(6);

ex = e*cos(w);
ey = e*sin(w);

eta = sqrt(1 - e^2);
kappa = 0.75 * J2*Re^2*sqrt(mu) / (a^3.5 * eta^4);
C = sin(w);
D = cos(w);
% E = 1+eta;
% F = 4+3*eta;
G = eta^-2;
% P = 3*cos(i)^2 - 1;
Q = 5*cos(i)^2 - 1;
% R = cos(i);
S = sin(2*i);
T = sin(i)^2;
% U = sin(i);
% V = tan(i/2);
% W = cos(i/2)^2;


A = kappa * ...
     [  0,          0,                  0,                  0,          0;     
        3.5*ey*Q,   -(4*ex*ey*G + C)*Q, -(1+4*ey^2*G-D)*Q,  5*ey*S,     0;  
        -3.5*ex*Q,  (1+4*ex^2*G-D)*Q,   (4*ex*ey*G - C)*Q,  -5*ex*S,    0;     
        0,          0,                  0,                  0,          0;
        3.5*S,      -4*ex*G*S,          -4*ey*G*S,          2*T,        0;
     ];

end