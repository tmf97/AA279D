function Phi = STM_J2(aoes, ti, tf, mu)
%STM_J2 Computes the State Transition matrix with J2 effects
%   For quasi-nonsingular ROEs. From Eqn 2.5 of Chernick
%   ti/ are the intial/final times in seconds.

% values hardcoded for Earth
J2 = 1.08263E-3; 
Re = 6378100; % m

a = aoes(1);
e = aoes(2);
i = aoes(3);
% Omega = aoes(4);
w = aoes(5);
% nu = aoes(6);

n = mean_motion(mu, a);
% eccentricity dependent parameters
ex = e*cos(w);
ey = e*sin(w);

eta = sqrt(1-e^2);
E = 1 + eta;
G = 1 / eta^2;
F = 4 + 3*eta;
gamma = (3/4) * J2*Re^2 * sqrt(mu);
kappa = gamma / (a^(7/2) * eta^4);


% inclination dependent parameters
P = 3*cos(i)^2 - 1;
Q = 5*cos(i)^2 - 1;
S = sin(2*i);
T = sin(i)^2;

tau = tf - ti;
wdot = kappa*Q;


% assuming exi = exf = ex, eyi = eyf = ey.

Phi = [
    1,  0,  0,  0,  0,  0;
    -7*kappa*eta*P*tau - (3/2)*n*tau, 1, 7*kappa*ex*P*tau/eta, 7*kappa*ey*P*tau/eta, -7*kappa*eta*S*tau, 0;
    3.5*kappa*ey*Q*tau,     0,  cos(wdot*tau) - 4*kappa*ex*ey*G*Q*tau, -sin(wdot*tau)-4*kappa*ey*ey*G*Q*tau, 5*kappa*ey*S*tau, 0;
    -3.5*kappa*ex*Q*tau,    0,  sin(wdot*tau) + 4*kappa*ex*ex*G*Q*tau, cos(wdot*tau)-4*kappa*ey*ex*G*Q*tau, -5*kappa*ex*S*tau, 0;
    0,  0,  0,  0,  1,  0;
    3.5*kappa*S*tau,    0, -4*kappa*ex*G*S*tau, -4*kappa*ey*G*S*tau, 2*kappa*T*tau, 1;
    ];
end