close all

% validate propagator:
MU = 3.9860043550702260E+14; % m^3/s^2
opts = odeset("RelTol",1e-13,"AbsTol",1e-15);


% load initial conditions
% mean AOE
maoe_chief_i = [10000000, 0.3, 0.9, 0.1, 0.1, 0.1];
a0 = maoe_chief_i(1);

mean_roes_initial = [0.1, 0.1, 100, 100, 0.1, 0.1].' / a0; 
mean_roes_desired = [0.1, 0.1, 10, 20, 0.1, 0.1].' / a0;

maoe_deputy_i = roe2aoe(maoe_chief_i, mean_roes_initial);


% Delta delta alpha
roe_difference = mean_roes_desired - mean_roes_initial;

% compute optimal maneuver locations
Ddey = roe_difference(4);
Ddex = roe_difference(3);
phi_ip = atan2(Ddey, Ddex);

Ddiy = roe_difference(6);
Ddix = roe_difference(5);
phi_oop = atan2(Ddiy, Ddix);

% compute delta-V lower bound
dex = roe_difference(3);
dey = roe_difference(4);
de = [dex; dey];

e = maoe_chief_i(2);
n = mean_motion(MU, a0);
eta = sqrt(1-e^2);

dv_LB = norm(de) * eta * n * a0 / sqrt(3*e^4 - 7*e^2 + 4);


% CONTROLLER GAINS
N = 12;
k = 5e2; % ?


dt = 1; % seconds
Tp = 2*pi/mean_motion(MU, a0);
time_range = 20*Tp;
num_steps = floor(time_range/dt);
tspan = linspace(0, time_range, num_steps);
dt = time_range / num_steps;

% initial conditions:
% convert mean to osculating
aoe_chief_i = mean2osc(maoe_chief_i, 1);
aoe_deputy_i = mean2osc(maoe_deputy_i, 1);

[rc0, vc0] = koe2pv(aoe_chief_i, MU);
[rd0, vd0] = koe2pv(aoe_deputy_i, MU);
initial_states_chief = [rc0; vc0];
initial_states_deputy = [rd0; vd0];

% setup state storage arrays
states_chief = zeros(6, length(tspan));
states_deputy = zeros(6, length(tspan));
koes_chief = zeros(6, length(tspan));
koes_deputy = zeros(6, length(tspan));
moes_chief = zeros(6, length(tspan));
moes_deputy = zeros(6, length(tspan));
mean_roes = zeros(6, length(tspan));
osc_roes = zeros(6, length(tspan));

delta_vs = zeros(3, length(tspan));

states_chief(:,1) = initial_states_chief;
states_deputy(:,1) = initial_states_deputy;
koes_chief(:,1) = aoe_chief_i;
koes_deputy(:,1) = aoe_deputy_i;
moes_chief(:,1) = maoe_chief_i;
moes_deputy(:,1) = maoe_deputy_i;
mean_roes(:,1) = mean_roes_initial;
osc_roes(:,1) = quasi_nonsingular_roe(aoe_chief_i, aoe_deputy_i);


for i=2:length(tspan)
    % compute Keplerian and J2 accelerations for both chief and deputy
    r_chief = states_chief(1:3, i-1);
    v_chief = states_chief(4:6, i-1);
    r_deputy = states_deputy(1:3, i-1);
    v_deputy = states_deputy(4:6, i-1);

    ds_J2_chief = newton3d_J2(0, [r_chief; v_chief]);
    ds_J2_deputy = newton3d_J2(0, [r_deputy; v_deputy]);

    a_chief = ds_J2_chief(4:6);
    a_deputy = ds_J2_deputy(4:6);

    % compute KOEs for both chief and deputy
    koe_chief = pv2koe([r_chief; v_chief], MU);
    koe_deputy = pv2koe([r_deputy; v_deputy], MU);

    moe_chief = osc2mean(koe_chief, 1);
    moe_deputy = osc2mean(koe_deputy, 1);
    roes = quasi_nonsingular_roe(moe_chief, moe_deputy);
    ddalpha = roes - mean_roes_desired;
   

    % store
    states_chief(:,i) = [r_chief; v_chief];
    states_deputy(:,i) = [r_deputy; v_deputy];
    koes_chief(:,i) = koe_chief;
    koes_deputy(:,i) = koe_deputy;
    moes_chief(:,i) = moe_chief;
    moes_deputy(:,i) = moe_deputy;
    mean_roes(:,i) = roes;

     % reform to delete dlambda
    roes = [roes(1); roes(3); roes(4); roes(5); roes(6)];
    ddalpha = [ddalpha(1); ddalpha(3); ddalpha(4); ddalpha(5); ddalpha(6)];


    e = koe_deputy(2);
    nu = koe_deputy(6);
    w = koe_deputy(5);
    M = true2mean(nu, e);

    % mean argument of latitude of the deputy
    phi = w + M;

    cosJ = cos(phi - phi_ip)^N;
    cosH = cos(phi - phi_oop)^N;

    % compute control law
    P = (1/k) * diag([cosJ, cosJ, cosJ, cosH, cosH]);
    A = plant_matrix_J2(koe_chief, MU);
    B = control_input_matrix_v2(koe_chief, MU);
    Binv = pinv(B);
    
    u_RTN = -Binv * (A * roes + P * ddalpha);
    u_RTN = [0; u_RTN(1); u_RTN(2)];

    % rotate control vector into ECI
    rhat = r_chief/norm(r_chief);
    h0 = cross(r_chief, v_chief);
    nhat = h0/norm(h0);
    that = cross(nhat, rhat);
    C_RTN_to_ECI = [rhat, that, nhat];

    u_ECI = C_RTN_to_ECI * u_RTN;

    delta_vs(:, i) = u_RTN;
    ds_chief = [v_chief; a_chief;];
    ds_deputy = [v_deputy; a_deputy + u_ECI];

    states_chief(:, i) = states_chief(:, i-1) + dt * ds_chief;
    states_deputy(:, i) = states_deputy(:, i-1) + dt * ds_deputy;
end

disp("Delta-V Used:");
disp(sum(vecnorm(delta_vs))*dt);

disp("Delta-V Lower-bound:");
disp(dv_LB);


% plot everything in RTN of the chief:
mms2_rtns = zeros(3, length(tspan));
for i=1:length(tspan)
    [mms2_rtns(:,i), ~] = pvs2rtn(states_chief(1:3,i),states_chief(4:6,i), states_deputy(1:3,i),states_deputy(4:6,i), MU, a0);
end

Rs2 = mms2_rtns(1, :);
Ts2 = mms2_rtns(2, :);
Ns2 = mms2_rtns(3, :);
tspan = tspan.';
figure
subplot(3,2,1);
hold on
plot(tspan / Tp, Rs2);
ylabel("R (m)")
grid on;

subplot(3,2,3);
hold on
plot(tspan / Tp, Ts2);
ylabel("T (m)")
grid on;


subplot(3,2,5);
hold on
plot(tspan / Tp, Ns2);
ylabel("N (m)")
xlabel("Orbital Periods")
grid on;


subplot(2,4,3)
hold on
plot(Ts2, Rs2)
xlabel("T (m)")
ylabel("R (m)")
grid on;
axis square;

subplot(2,4,4)
hold on
plot(Ns2, Rs2)
xlabel("N (m)")
ylabel("R (m)")
grid on;
axis square;

subplot(2,4,7)
hold on
plot(Ts2, Ns2)
xlabel("T (m)")
ylabel("N (m)")
grid on;
axis square;

subplot(2,4,8)
hold on
plot3(Rs2, Ts2, Ns2)
xlabel("R (m)")
ylabel("T (m)")
zlabel("N (m)")
grid on;
axis square;

% plot deltaV
figure
subplot(2,1,1)
hold on 
plot(tspan/Tp, dt*delta_vs(1,:), 'DisplayName','R')
plot(tspan/Tp, dt*delta_vs(2,:), 'DisplayName','T')
plot(tspan/Tp, dt*delta_vs(3,:), 'DisplayName','N')
ylabel('Delta V (m/s)')
hold off
legend
grid on

subplot(2,1,2)
plot(tspan/Tp, cumsum(dt*vecnorm(delta_vs)), 'DisplayName','Total DeltaV')
xlabel('Orbit Number')
ylabel('Delta V (m/s)')
grid on

% plot ROEs
figure
subplot(1,3,1)
hold on
plot(a0*mean_roes(2, :), a0*mean_roes(1, :), 'DisplayName', 'Trajectory')
scatter(a0*mean_roes_desired(2), a0*mean_roes_desired(1), 'DisplayName', 'Target')
scatter(a0*mean_roes(2, 1), a0*mean_roes(1, 1), 'DisplayName','Start');
xlabel("a\delta\lambda (m)")
ylabel("a\delta a (m)")
grid on;
axis square;

hold off
legend


subplot(1,3,2)
hold on
plot(a0*mean_roes(3, :), a0*mean_roes(4, :))
scatter(a0*mean_roes_desired(3), a0*mean_roes_desired(4))
scatter(a0*mean_roes(3, 1), a0*mean_roes(4, 1));
xlabel("a\delta e_x (m)")
ylabel("a\delta e_y (m)")
grid on;
axis square;

hold off

subplot(1,3,3)
hold on
plot(a0*mean_roes(5, :), a0*mean_roes(6, :))
scatter(a0*mean_roes_desired(5), a0*mean_roes_desired(6))
scatter(a0*mean_roes(5, 1), a0*mean_roes(6, 1));
xlabel("a\delta i_x (m)")
ylabel("a\delta i_y (m)")
grid on;
axis square;

hold off

