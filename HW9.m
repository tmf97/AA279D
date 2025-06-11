% HW9
% Oh yeah, it's all coming together
% Constants
J2 = 0.0010826358191967; % Earth's J2 coefficient
mu = 3.986004415e14; % Earth gravitational parameter [m^3/s^2]
rE = 6.378136300e6; % Earth radius [m]

% Plotting
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesFontSize',18);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Problem 2b
% Alter the control law state input by corrupting the ground truth with some
% representative noise. This is the same process you used to generate your first
% set of measurements in Problem Set 7-8. How does your controller perform?
% What are the main differences compared to the original implementation?
% References: Chernick PhD
% Starter code is NOT self-encompassing. You are expected to implement some
% additional functions. Everything is `impulsive_control.m` should have
% function headers at a minimum though. 

close all; clear; clc;

% Integration
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial chief absolute singular orbital parameters
sma_c  = 68920.927e3;      % semi-major axis [m] %sma_c = 10000e3;
ecc_c  = 1e-3;            % eccentricity component in I
inc_c  = deg2rad(10);  % inclination [rad]
raan_c = deg2rad(10);    % RAAN [rad]
aop_c  = deg2rad(10);      % aop [rad]
M_c    = deg2rad(0);      % mean anomaly [rad]

oe_init_c = [sma_c, ecc_c, inc_c, raan_c, aop_c, M_c].'; % combine the above into a single vector
% oe_init_c_j2 = mean2osc(oe_init_c, 1);

[r_init_c,v_init_c] = koe2pv(oe_init_c, mu); % chief position and velocity in inertial frame, J2-perturbed
rv_init_c = [r_init_c;v_init_c];

% Helpful parameters
T = 2*pi*sqrt(sma_c^3/mu); % [sec]
n = sqrt(mu/sma_c^3); % mean motion

% Initial conditions of the Deputy as QNS ROEs [dsma;dlambda;dex;dey;dix;diy]
aroe_init1 = [0, 0, 50, 100,0,0].'; % m
roe_init1 = aroe_init1 / sma_c;

% Initial absolute conditions of the deputy
oe_init_d1 = roe2aoe(oe_init_c,roe_init1);

[r_init_d,v_init_d] = koe2pv(oe_init_d1, mu);
rv_init_d = [r_init_d;v_init_d];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconfiguration Conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desired ROE state 
aroe_des = [1000, 0, 50, 100, 0, 0].'; % m
% aroe_des = [0, 0, 50, 1000, 0, 0].'; % m

% Desired reconfiguration time
dt = 2*T; % 10 orbit periods

% Simulation step size
sim_step = 100; % s

%% Corrupt initial state estimate
% but not the true dynamics
aroe_init1 = aroe_init1 + [100,100,100,100,100,50].';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine optimal maneuver times and magnitudes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the functions you'd like to use for control input matrix and STM
stm = @(chief_oe, t) chernick_J2_stm(chief_oe, t, rE, mu, J2); 
control_input_matrix = @(chief_oe) chernick_control_matrix(chief_oe, mu);

% Get the optimal maneuver plan for the reconfiguration window
[t_maneuvers, maneuvers, total_cost] = impulsive_control(oe_init_c, aroe_init1, aroe_des, dt, stm, control_input_matrix, rE, mu, J2);
aDdalpha = aroe_des - aroe_init1;
dv_LB = dv_lower_bound(aDdalpha,sma_c,ecc_c, dt,mu);

disp("Total Delta-V Cost (m/s):")
disp(total_cost)
disp("Delta-V Lower Bound(m/s):")
disp(dv_LB)

% Round the time
t_maneuvers = round(t_maneuvers./sim_step).*sim_step;

%% Simulate!
% Finally!
opts = odeset("RelTol",1e-13,"AbsTol",1e-15);

states_chief = rv_init_c;
states_deputy = rv_init_d;

propagation_times_start = [0, t_maneuvers];
propagation_times_end = [t_maneuvers, dt];
times = [];
dvs = [];
for i = 1:length(t_maneuvers)+1
    t_start = propagation_times_start(i);
    t_end = propagation_times_end(i);
    ts = t_start:sim_step:t_end;
    assert(t_end > t_start)
    
    % propagate chief
    s0c = states_chief(:, end);
    [prop_ts, propagated_states] = ode45(@newton3d_J2, ts, s0c, opts);
    states_chief = [states_chief, propagated_states.'];
    times = [times; prop_ts];
    % propagate deputy
    s0d = states_deputy(:, end);
    [~, propagated_states] = ode45(@newton3d_J2, ts, s0d, opts);
    states_deputy = [states_deputy, propagated_states.'];

    % apply maneuver dv
    if i <= length(t_maneuvers)    
        last_state = states_deputy(:, end);
        r_last = last_state(1:3);
        v_last = last_state(4:6);
        dv_RTN = maneuvers(:, i);
        % rotate vector:
        dv_ECI = rtn2eci(r_last, v_last, dv_RTN);
        state_override = [0; 0; 0; dv_ECI];
        states_deputy(:, end) = states_deputy(:, end) + state_override;
        states_chief(:, end) = states_chief(:, end);
    else
        states_deputy(:, end) = [];
        states_chief(:, end) = [];
    end

    % dv accounting:
    if i > 1
        prior_dv = dvs(end);
        dv_mag = norm(dv_RTN);
    else
        dv_mag = 0;
        prior_dv = 0;
    end

    dvs = [dvs, (prior_dv+dv_mag)*ones(1,length(prop_ts))];
end


% plot everything in RTN of the chief:
mms2_rtns = zeros(3, length(states_chief));
mms2_rtn_vs = zeros(3, length(states_chief));
mean_roes = zeros(6, length(states_chief));
for i=1:length(states_deputy)
    r_chief = states_chief(1:3,i);
    v_chief = states_chief(4:6,i);

    r_deputy = states_deputy(1:3,i);
    v_deputy = states_deputy(4:6,i);
    [mms2_rtns(:,i), mms2_rtn_vs(:,i)] = pvs2rtn(r_chief, v_chief, r_deputy, v_deputy);
    % compute KOEs for both chief and deputy
    koe_chief = pv2koe([r_chief; v_chief], mu);
    koe_deputy = pv2koe([r_deputy; v_deputy], mu);

    moe_chief = osc2mean(koe_chief, 1);
    moe_deputy = osc2mean(koe_deputy, 1);
    mean_roes(:, i) = quasi_nonsingular_roe(moe_chief, moe_deputy);
end

Rs2 = mms2_rtns(1, :);
Ts2 = mms2_rtns(2, :);
Ns2 = mms2_rtns(3, :);
times = times.';
times = times / T;

figure
subplot(3,2,1);
hold on
plot(times, Rs2);
ylabel("R (m)")
grid on;

subplot(3,2,3);
hold on
plot(times, Ts2);
ylabel("T (m)")
grid on;


subplot(3,2,5);
hold on
plot(times, Ns2);
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
axis square equal;
pbaspect([1,1,1])



% % 3D!
% figure
% grid on
% hold on
% axis equal
% Re = 6378000; % m
% [X,Y,Z] = sphere;
% surf(X*Re, Y*Re, Z*Re, 'DisplayName', 'Earth')
% plot3(states_chief(1,:), states_chief(2,:), states_chief(3,:))
% plot3(states_deputy(1,:), states_deputy(2,:), states_deputy(3,:))


% plot ROEs
figure
subplot(1,3,1)
hold on
plot(sma_c*mean_roes(2, :), sma_c*mean_roes(1, :), 'DisplayName', 'Trajectory')
scatter(aroe_des(2), aroe_des(1), 'DisplayName', 'Target')
scatter(sma_c*mean_roes(2, 1), sma_c*mean_roes(1, 1), 'DisplayName','Start');
xlabel("$a \delta \lambda$ (m)")
ylabel("$a \delta a$ (m)")
grid on;
axis square;

hold off
legend


subplot(1,3,2)
hold on
plot(sma_c*mean_roes(3, :), sma_c*mean_roes(4, :))
scatter(aroe_des(3), aroe_des(4))
scatter(sma_c*mean_roes(3, 1), sma_c*mean_roes(4, 1));
xlabel("$a \delta e_x$ (m)")
ylabel("$a \delta e_y$ (m)")
grid on;
axis square;

hold off

subplot(1,3,3)
hold on
plot(sma_c*mean_roes(5, :), sma_c*mean_roes(6, :))
scatter(aroe_des(5), aroe_des(6))
scatter(sma_c*mean_roes(5, 1), sma_c*mean_roes(6, 1));
xlabel("$a \delta i_x$ (m)")
ylabel("$a \delta i_y$ (m)")
grid on;
axis square;

hold off

figure
subplot(2,1,1)
hold on
plot(times, abs(sma_c*mean_roes(1, :) - aroe_des(1)), "DisplayName", "$\delta a$")
plot(times, abs(sma_c*mean_roes(2, :) - aroe_des(2)), "DisplayName", "$\delta\lambda$")
plot(times, abs(sma_c*mean_roes(3, :) - aroe_des(3)), "DisplayName", "$\delta e_x$")
plot(times, abs(sma_c*mean_roes(4, :) - aroe_des(4)), "DisplayName", "$\delta e_y$")
plot(times, abs(sma_c*mean_roes(5, :) - aroe_des(5)), "DisplayName", "$\delta i_x$")
plot(times, abs(sma_c*mean_roes(6, :) - aroe_des(6)), "DisplayName", "$\delta i_y$")
legend('Location','eastoutside')
xlabel("Orbital Periods")
ylabel("ROE Error Magnitude (m)")
ylim([1e-5, 1e5]);
yscale log
grid on;
hold off

subplot(2,1,2);
hold on
plot(times, dvs, "DisplayName","Controller")
yline(dv_LB, "DisplayName","Lower Bound")
xlabel("Orbital Periods")
ylabel("$\Delta V$ Consumed (m/s)")
legend('Location','eastoutside')
grid on
hold off


