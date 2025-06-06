function [Phi] = chernick_J2_stm(singular_oe, t, rP, mu, J2)
%{
    chernick_J2_stm calculates the J2 STM for ROEs defined in Chernick 2.4 
    
    Inputs:
        - singular_oe: keplerian orbit elements [a, e, i, RAAN, w, M] 
        - t: time at which to evaluate the STM (s)
        - Re: radius of the primary attractor (m)
        - mu: Gravitational parameter of the primary attractor (m^3/s^2)
        - J2: J2 coefficient of the primary attractor
        
    Outputs: 
        - Phi: STM evaluated at t
%}

a = singular_oe(1);
e = singular_oe(2);
i = singular_oe(3);
% Omega = singular_oe(4);
w = singular_oe(5);
% nu = singular_oe(6);

n = mean_motion(mu, a);
% eccentricity dependent parameters
ex = e*cos(w);
ey = e*sin(w);

eta = sqrt(1-e^2);
E = 1 + eta;
G = 1 / eta^2;
F = 4 + 3*eta;
gamma = (3/4) * J2*rP^2 * sqrt(mu);
kappa = gamma / (a^(7/2) * eta^4);


% inclination dependent parameters
P = 3*cos(i)^2 - 1;
Q = 5*cos(i)^2 - 1;
S = sin(2*i);
T = sin(i)^2;

wdot = kappa*Q;


% assuming exi = exf = ex, eyi = eyf = ey.

Phi = [
    1,                              0,  0,                                  0,                                  0,                  0;
    -7*kappa*eta*P*t - (3/2)*n*t,   1,  7*kappa*ex*P*t/eta,                 7*kappa*ey*P*t/eta,                 -7*kappa*eta*S*t,   0;
    3.5*kappa*ey*Q*t,               0,  cos(wdot*t) - 4*kappa*ex*ey*G*Q*t,  -sin(wdot*t)-4*kappa*ey*ey*G*Q*t,   5*kappa*ey*S*t,     0;
    -3.5*kappa*ex*Q*t,              0,  sin(wdot*t) + 4*kappa*ex*ex*G*Q*t,  cos(wdot*t)-4*kappa*ey*ex*G*Q*t,    -5*kappa*ex*S*t,    0;
    0,                              0,  0,                                  0,                                  1,                  0;
    3.5*kappa*S*t,                  0,  -4*kappa*ex*G*S*t,                  -4*kappa*ey*G*S*t,                  2*kappa*T*t,        1;
    ];

% assume near-circular orbits for now:
% Phi = [ 1,                      0,  0,              0,              0,              0;
%         -(7*kappa*E*P+3*n)*t/2, 1,  0,              0,              0,              0;
%         0,                      0,  cos(wdot*t),    -sin(wdot*t),   0,              0;
%         0,                      0,  sin(wdot*t),    cos(wdot*t),    0,              0;
%         0,                      0,  0,              0,              1,              0;
%         3.5*kappa*S*t,          0,  0,              0,              2*kappa*T*t,    1;
% ];
end

