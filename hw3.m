close all
MU = 3.9860043550702260E+14; % m^3/s^2

% initial Keplerian Orbital elements of Chief spacecraft (FMMS-1):
a0 = 63781000; % semi-major axis, meters
e0 = 1e-5; % eccentricity
i0 = 6.311884598781460E+01; % inclination, degrees
raan0 = 3.567045162674137E+02; % Right Ascencion of the Ascending Node (aka longitude of ascending node aka omega), degrees
argp0 = 1.582342516847609E+02; % Argument of periapsis (aka Argument of Perifocus), degrees
nu0 = 0; % True Anomaly, degrees

i0 = deg2rad(i0);
raan0 = deg2rad(raan0);
argp0 = deg2rad(argp0);
nu0 = deg2rad(nu0);

% initial Keplerian Orbital elements of deputy spacecraft (FMMS-2)
a1 = a0; % semi-major axis
e1 = e0; % eccentricity
i1 =  i0 - deg2rad(0.01); % inclination
raan1 = raan0; %Right Ascencion of the Ascending Node
argp1 = argp0; % Argument of periapsis
nu1 = deg2rad(0.01); % True Anomaly, degrees

% check that the difference is less than 0.001 of r0
[r0, v0] = koe2pv(a0, e0, rad2deg(i0), rad2deg(raan0), rad2deg(argp0), rad2deg(nu0));
[r1, v1] = koe2pv(a1, e1, rad2deg(i1), rad2deg(raan1), rad2deg(argp1), rad2deg(nu1));


rho_ECI = r1-r0;
disp("rho/r0:");
disp(norm(rho_ECI/norm(r0)));
assert(norm(rho_ECI/norm(r0)) < 1e-3, "Rho too high!")

disp("Initial Position and Velocities in ECI");
disp("    Sat Name    x (m)               y (m)             z (m)              vx (m/s)       vy (m/s)        vz (m/s)")
disp(["FMMS-1", r0.', v0.'])
disp(["FMMS-2", r1.', v1.'])


% define initial RTN frame
rhat = r0/norm(r0);
h0 = cross(r0, v0);
nhat = h0/norm(h0);
that = cross(nhat, rhat);
C_ECI_to_RTN = [rhat, that, nhat].';

% angular speed of the RTN frame in RTN

% fdot is the time derivative of true anomaly.
% r^2 * fdot = h
fdot = norm(h0) / (norm(r0)^2);

% sanity check this against mean motion, should be similar but not quite
% idential
n = mean_motion(MU, a0);
assert(abs((fdot-n)/n) < 1e-3, "Are you sure fdot is correct?");
omega_RTN_in_ECI_in_RTN = [0;0; fdot];

rho_RTN = C_ECI_to_RTN * rho_ECI;

rho_dot_ECI = v1 - v0;

omega_ECI_in_RTN_in_ECI = C_ECI_to_RTN.' * -omega_RTN_in_ECI_in_RTN;

rho_dot_RTN_in_ECI = rho_dot_ECI + cross(omega_ECI_in_RTN_in_ECI, rho_ECI);
rho_dot_RTN_in_RTN = C_ECI_to_RTN * rho_dot_RTN_in_ECI;

disp("rho RTN (m):")
disp(rho_RTN)

disp("rho dot RTN (m/s):")
disp(rho_dot_RTN_in_RTN)


% compute HCW integration constants
x0 = [rho_RTN; rho_dot_RTN_in_RTN];
an = sqrt(MU/a0);
A0 = diag([a0, a0, a0, an, an, an]);
% assuming t = 0
function A = HCW_matrix(t, n)
A = [1,         sin(n*t),   cos(n*t),       0, 0,           0;
    -(3/2)*n*t, 2*cos(n*t), -2*sin(n*t),    1, 0,           0;
    0,          0,          0,              0, sin(n*t),    cos(n*t);
    0,          cos(n*t),   -sin(n*t),      0, 0,           0;
    -3/2,       -2*sin(n*t),-2*cos(n*t),    0, 0,           0;
    0,          0,          0,              0, cos(n*t),    -sin(n*t)];
end

K = (A0 * HCW_matrix(0, mean_motion(MU, a0))) \ x0;

disp("K:")
disp(K)


% Sanity check:
% Page 86 if Spacecraft Formation Flying
function states = solve_HCW(initial_rho, initial_rho_dot, mean_motion, t)
x = initial_rho(1);
y = initial_rho(2);
z = initial_rho(3);
xdot = initial_rho_dot(1);
ydot = initial_rho_dot(2);
zdot = initial_rho_dot(3);
n = mean_motion;

xs = (4*x + 2 * ydot / n) + (xdot/n) * sin(n*t) - (3*x + 2*ydot/n) * cos(n*t);
ys = -(6*n*x + 3*ydot)*t + (y - 2*xdot/n) + (6*x + 4*ydot/n)*sin(n*t) + 2*xdot * cos(n*t)/n;
zs = zdot * sin(n*t)/n + z*cos(n*t);

states = [xs; ys; zs];
end

Tp = 2*pi*sqrt(a0^3 / MU);
t = linspace(0, 15*Tp, 1000);
states = solve_HCW(rho_RTN, rho_dot_RTN_in_RTN, n, t);

Rs_th = states(1, :);
Ts_th = states(2, :);
Ns_th = states(3, :);

figure
plot_t = t / Tp;

subplot(3,2,1);
plot(plot_t, Rs_th);
ylabel("R (m)")
grid on;


subplot(3,2,3);
plot(plot_t, Ts_th);
ylabel("T (m)")
grid on;


subplot(3,2,5);
plot(plot_t, Ns_th);
ylabel("N (m)")
xlabel("Orbital Periods")
grid on;


subplot(2,4,3)
plot(Ts_th, Rs_th)
xlabel("T (m)")
ylabel("R (m)")
grid on;
axis square;

subplot(2,4,4)
plot(Ns_th, Rs_th)
xlabel("N (m)")
ylabel("R (m)")
grid on;
axis square;

subplot(2,4,7)
plot(Ts_th, Ns_th)
xlabel("T (m)")
ylabel("N (m)")
grid on;
axis square;

subplot(2,4,8)
plot3(Rs_th, Ts_th, Ns_th)
xlabel("R (m)")
ylabel("T (m)")
zlabel("N (m)")
grid on;
axis square;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   PROBLEM 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Redefine intial eccentricity
e0 = 0.1;
e1 = e0;

% check that the difference is less than 0.001 of r0
[r0, v0] = koe2pv(a0, e0, rad2deg(i0), rad2deg(raan0), rad2deg(argp0), rad2deg(nu0));
[r1, v1] = koe2pv(a1, e1, rad2deg(i1), rad2deg(raan1), rad2deg(argp1), rad2deg(nu1));


rho_ECI = r1-r0;
disp("rho/r0:");
disp(norm(rho_ECI/norm(r0)));
assert(norm(rho_ECI/norm(r0)) < 1e-3, "Rho too high!")

disp("Initial Position and Velocities in ECI");
disp("    Sat Name    x (m)               y (m)             z (m)              vx (m/s)       vy (m/s)        vz (m/s)")
disp(["FMMS-1", r0.', v0.'])
disp(["FMMS-2", r1.', v1.'])


% define initial RTN frame
rhat = r0/norm(r0);
h0 = cross(r0, v0);
nhat = h0/norm(h0);
that = cross(nhat, rhat);
C_ECI_to_RTN = [rhat, that, nhat].';

% angular speed of the RTN frame in RTN

% fdot is the time derivative of true anomaly.
% r^2 * fdot = h
fdot = norm(h0) / (norm(r0)^2);

% sanity check this against mean motion, should be similar but not quite
% idential
n = mean_motion(MU, a0);
assert(abs((fdot-n)/n) < 0.5, "Are you sure fdot is correct?");
omega_RTN_in_ECI_in_RTN = [0;0; fdot];

rho_RTN = C_ECI_to_RTN * rho_ECI;

rho_dot_ECI = v1 - v0;

omega_ECI_in_RTN_in_ECI = C_ECI_to_RTN.' * -omega_RTN_in_ECI_in_RTN;

rho_dot_RTN_in_ECI = rho_dot_ECI + cross(omega_ECI_in_RTN_in_ECI, rho_ECI);
rho_dot_RTN_in_RTN = C_ECI_to_RTN * rho_dot_RTN_in_ECI;

disp("rho RTN (m):")
disp(rho_RTN)

disp("rho dot RTN (m/s):")
disp(rho_dot_RTN_in_RTN)

% now normalize initial conditions
rho_RTN_norm = rho_RTN / norm(r0);
rho_dot_RTN_norm = rho_dot_RTN_in_RTN * fdot;

assert(norm(rho_RTN_norm) < 1e-3, "Rho too high!")


function states = tschauner_hempel(integration_constants, true_anoms, eccentricity, sma, h)
MU = 3.9860043550702260E+14; % m^3/s^2, for Earth

c1 = integration_constants(1);
c2 = integration_constants(2);
c3 = integration_constants(3);
c4 = integration_constants(4);
c5 = integration_constants(5);
c6 = integration_constants(6);
e = eccentricity;
fs = true_anoms;
weird = fs - mod(fs, 2*pi);
Ms = true2mean(fs, e);
Ms = Ms + weird;
n = mean_motion(MU, sma);
ts = Ms / n;
ts = ts - ts(1);

ks = 1+e*cos(fs);
Is = MU^2 * ts / norm(h)^3;
xbars = c1 * ks .* sin(fs) + c2*ks.*cos(fs) + c3*(2-3*e*ks.*Is.*sin(fs));
ybars = c4 + c1*(1+1./ks).*cos(fs) -c2*(1+1./ks).*sin(fs) - 3*c3*ks.^2.*Is;
zbars = c5 * cos(fs) + c6*sin(fs);

states = [xbars; ybars; zbars; ts];
end

Tp = 2*pi*sqrt(a0^3 / MU);
fs = linspace(0, 15*2*pi, 1000) + nu0;

initial_conditions = [rho_RTN_norm; rho_dot_RTN_norm];

% YA matrix for t=0
k0 = 1 + e0 * cos(nu0);
ksinf_prime = cos(nu0) + e0*cos(2*nu0);
kcosf_prime = -(sin(nu0) + e0*sin(2*nu0));

mat = [1,   k0*sin(nu0),    k0*cos(nu0), 0, 0, 0;
    0,   (1+k0)*cos(nu0), -(1+k0)*sin(nu0), 1, 0, 0;
    0,   0,              0,      0,      sin(nu0), cos(nu0);
    0,   ksinf_prime,    kcosf_prime,    0, 0, 0;
    -1.5, -2*k0*sin(nu0), e0 - 2*k0*cos(nu0), 0, 0, 0;
    0,   0,              0,      0,      cos(nu0), -sin(nu0)];


Cs = mat \ initial_conditions;
disp("YA integration Constants:");
disp(Cs);
states_th = tschauner_hempel(Cs,fs,e0,a0,h0);

Rs_th = states_th(1, :);
Ts_th = states_th(2, :);
Ns_th = states_th(3, :);
times = states_th(4, :);


figure
subplot(3,2,1);
plot(fs, Rs_th);
ylabel("R")
grid on;


subplot(3,2,3);
plot(fs, Ts_th);
ylabel("T")
grid on;


subplot(3,2,5);
plot(fs, Ns_th);
ylabel("N")
xlabel("True Anomaly [rad]")
grid on;


subplot(2,4,3)
plot(Ts_th, Rs_th)
xlabel("T")
ylabel("R")
grid on;
axis square;

subplot(2,4,4)
plot(Ns_th, Rs_th)
xlabel("N")
ylabel("R")
grid on;
axis square;

subplot(2,4,7)
plot(Ts_th, Ns_th)
xlabel("T")
ylabel("N")
grid on;
axis square;

subplot(2,4,8)
plot3(Rs_th, Ts_th, Ns_th)
xlabel("R")
ylabel("T")
zlabel("N")
grid on;
axis square;

% the trends in R and T are expected from the solution of the TH Equations:
% The I term grows with time, and so the magnitude of the oscillations of R
% will increase for a given non-zero c3 value given an eccentric orbit.
% Same goes for the average value of T,
% because its average value changes linearly with I.
% The relative motion in R and T is not bounded, which does not match our
% expectation from da = 0, which would presume a constant rho over each
% orbit. The motion is unbounded because the c3 integration constant is
% non-zero. We can bound the motion by choosing a set of initial conditions
% to set c3 to zero. The constraints on the initial conditions required to
% get c3 = 0, comes from k(f(0))y'(f(0)) + esin(f(0))[x'(f(0))-y(f(0))] +
% [2+ecosf(0)]x(f(0)) = 0, where k = 1 + e cosf. Note that these quantities
% are not normalized.


% part e - quasi-nonsingular ROE


disp("Quasi-nonsingular ROEs:")
M0 = true2mean(nu0, e0);
M1 = true2mean(nu1, e1);
QNS_ROEs = quasi_nonsingular_roe(a0, M0, argp0, e0, raan0, i0, a1, M1, argp1, e1, raan1, i1);
disp(QNS_ROEs)
% only the dlambda and dix elements are non-zero, because the two orbits
% only have inclination and true anomaly offsets


function states = linear_mapping(ROEs, nus, argp, e, i)
da = ROEs(1);
dlambda = ROEs(2);
dex = ROEs(3);
dey = ROEs(4);
dix = ROEs(5);
diy = ROEs(6);

eta = sqrt(1-e^2);
ex = e*cos(argp);
ey = e*sin(argp);
us = nus + argp;
ks = 1 + ex * cos(us) + ey * sin(us);
kprimes = -ex * sin(us) + ey * cos(us);

x_bar = da + (ks / eta^3) .* ( ...
    - (kprimes*dlambda) ...
    - (dex*cos(us) - dey*sin(us)) ...
    + ((ks-1)./(1+eta)) * (ex*dex + ey*dey) ...
    + kprimes*diy*cot(i));

y_bar = ((ks.^2 * dlambda / eta^3) ...
    + dex*(1+ks).*sin(us)/eta^2 ...
    - dey*(1+ks).*cos(us)/eta^2 ...
    + (1/eta^3)*(eta + ks.^2/(1+eta))*(ey*dex - ex*dey) ...
    + (1-ks.^2/eta^3)*diy*cot(i));

z_bar = dix*sin(us) - diy*cos(us);

states = [x_bar; y_bar; z_bar;];
end

states_lm = linear_mapping(QNS_ROEs, fs, argp0, e0, i0);

Rs_lm = states_lm(1, :);
Ts_lm = states_lm(2, :);
Ns_lm = states_lm(3, :);


figure
subplot(3,2,1);
hold on
plot(fs, Rs_th, 'DisplayName', 'TH');
plot(fs, Rs_lm, 'DisplayName', 'LM');
ylabel("R")
grid on;
hold off


subplot(3,2,3);
hold on
plot(fs, Ts_th, 'DisplayName', 'TH');
plot(fs, Ts_lm, 'DisplayName', 'LM');
ylabel("T")
grid on;


subplot(3,2,5);
hold on
plot(fs, Ns_th, 'DisplayName', 'TH');
plot(fs, Ns_lm, 'DisplayName', 'LM');
ylabel("N")
xlabel("True Anomaly [rad]")
grid on;


subplot(2,4,3)
hold on
plot(Ts_th, Rs_th, 'DisplayName', 'TH');
plot(Ts_lm, Rs_lm, 'DisplayName', 'LM');
xlabel("T")
ylabel("R")
grid on;
axis square;

subplot(2,4,4)
hold on
plot(Ns_th, Rs_th, 'DisplayName', 'TH');
plot(Ns_lm, Rs_lm, 'DisplayName', 'LM');
xlabel("N")
ylabel("R")
grid on;
axis square;

subplot(2,4,7)
hold on
plot(Ts_th, Ns_th, 'DisplayName', 'TH');
plot(Ts_lm, Ns_lm, 'DisplayName', 'LM');
xlabel("T")
ylabel("N")
grid on;
axis square;

subplot(2,4,8)
hold on
plot3(Rs_th, Ts_th, Ns_th, 'DisplayName', 'TH');
plot3(Rs_lm, Ts_lm, Ns_lm, 'DisplayName', 'LM');
xlabel("R")
ylabel("T")
zlabel("N")
grid on;
axis square;
legend



% do it again but keplerian
nus0 = keplerian_propagation(nu0, e0, a0, times);
nus1 = keplerian_propagation(nu1, e1, a1, times);

Rs_kep = zeros(size(times));
Ts_kep = zeros(size(times));
Ns_kep = zeros(size(times));
% spaghetti copy/paste code from above because I'm cooked
for i = 1:length(t)
[r0, v0] = koe2pv(a0, e0, rad2deg(i0), rad2deg(raan0), rad2deg(argp0), rad2deg(nus0(i)));
[r1, v1] = koe2pv(a1, e1, rad2deg(i1), rad2deg(raan1), rad2deg(argp1), rad2deg(nus1(i)));

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
% sanity check this against mean motion, should be similar but not quite
% idential
n = mean_motion(MU, a0);
assert(abs((fdot-n)/n) < 0.5, "Are you sure fdot is correct?");
omega_RTN_in_ECI_in_RTN = [0;0; fdot];

rho_RTN = C_ECI_to_RTN * rho_ECI;

rho_dot_ECI = v1 - v0;

omega_ECI_in_RTN_in_ECI = C_ECI_to_RTN.' * -omega_RTN_in_ECI_in_RTN;

rho_dot_RTN_in_ECI = rho_dot_ECI + cross(omega_ECI_in_RTN_in_ECI, rho_ECI);
rho_dot_RTN_in_RTN = C_ECI_to_RTN * rho_dot_RTN_in_ECI;

rho_RTN_normalized = rho_RTN / norm(r0);
Rs_kep(i) = rho_RTN_normalized(1);
Ts_kep(i) = rho_RTN_normalized(2);
Ns_kep(i) = rho_RTN_normalized(3);
end


figure
subplot(3,2,1);
hold on
plot(fs, Rs_th, 'DisplayName', 'TH');
plot(fs, Rs_lm, 'DisplayName', 'LM');
plot(fs, Rs_kep, 'DisplayName', 'Keplerian');

ylabel("R")
grid on;
hold off


subplot(3,2,3);
hold on
plot(fs, Ts_th, 'DisplayName', 'TH');
plot(fs, Ts_lm, 'DisplayName', 'LM');
plot(fs, Ts_kep, 'DisplayName', 'Keplerian');
ylabel("T")
grid on;


subplot(3,2,5);
hold on
plot(fs, Ns_th, 'DisplayName', 'TH');
plot(fs, Ns_lm, 'DisplayName', 'LM');
plot(fs, Ns_kep, 'DisplayName', 'Keplerian');

ylabel("N")
xlabel("True Anomaly [rad]")
grid on;


subplot(2,4,3)
hold on
plot(Ts_th, Rs_th, 'DisplayName', 'TH');
plot(Ts_lm, Rs_lm, 'DisplayName', 'LM');
plot(Ts_kep, Rs_kep, 'DisplayName', 'Keplerian');
xlabel("T")
ylabel("R")
grid on;
axis square;

subplot(2,4,4)
hold on
plot(Ns_th, Rs_th, 'DisplayName', 'TH');
plot(Ns_lm, Rs_lm, 'DisplayName', 'LM');
plot(Ns_kep, Rs_kep, 'DisplayName', 'Keplerian');
xlabel("N")
ylabel("R")
grid on;
axis square;

subplot(2,4,7)
hold on
plot(Ts_th, Ns_th, 'DisplayName', 'TH');
plot(Ts_lm, Ns_lm, 'DisplayName', 'LM');
plot(Ts_kep, Ns_kep, 'DisplayName', 'Keplerian');
xlabel("T")
ylabel("N")
grid on;
axis square;

subplot(2,4,8)
hold on
plot3(Rs_th, Ts_th, Ns_th, 'DisplayName', 'TH');
plot3(Rs_lm, Ts_lm, Ns_lm, 'DisplayName', 'LM');
plot3(Rs_kep, Ts_kep, Ns_kep, 'DisplayName', 'Keplerian');
xlabel("R")
ylabel("T")
zlabel("N")
grid on;
axis square;
legend


% plot errors
figure
subplot(2,1,1)
x = fs / (2 * pi);
hold on 
plot(x, Rs_lm - Rs_kep, 'DisplayName', 'R')
plot(x, Ts_lm - Ts_kep, 'DisplayName', 'T')
plot(x, Ns_lm - Ns_kep, 'DisplayName', 'N')
ylabel("Error (GLM)")



subplot(2,1,2)
hold on 

legend
ylabel("Error (T-H)")
plot(x, Rs_th - Rs_kep, 'DisplayName', 'R')
plot(x, Ts_th - Ts_kep, 'DisplayName', 'T')
plot(x, Ns_th - Ns_kep, 'DisplayName', 'N')
xlabel('Orbit Number')