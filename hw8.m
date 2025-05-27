close all
% Kalman Filter time, baby!
% constants
MU = 3.9860043550702260E+14; % m^3/s^2
opts = odeset("RelTol",1e-11,"AbsTol",1e-13);


% initial absolute orbital elements
a0 = 10000000;
e = 0.85;
inc = deg2rad(28.5);
Omega = deg2rad(15);
w = deg2rad(15);
nu = deg2rad(90);
[r0, v0] = koe2pv([a0, e, inc, Omega, w, nu], MU);
initial_states_chief = [r0; v0];


Tp = 2*pi / mean_motion(MU, a0);
num_steps = 1e4;
ts = linspace(0, 2*Tp, num_steps);

% setup state storage arrays
truth_states_chief = zeros(6, length(ts));
pre_measurement_state_estimate = zeros(6, length(ts));
truth_states_chief(:,1) = initial_states_chief;

estimated_states_chief = zeros(6, length(ts));

initial_offsets = [50000; -5000; 5000; 10; 10; 10];
measurement_noise = [500; 500; 500; 0.1; 0.1; 0.1];
estimated_states_chief(:,1) = initial_states_chief + initial_offsets;
pre_measurement_state_estimate(:,1) = initial_states_chief + initial_offsets;


Ps = zeros(6,6,length(ts));
P0 = diag(abs(initial_offsets));
Ps(:,:,1) = P0;

for i=2:length(ts)
    dt_now = ts(i) - ts(i-1);
    % get true anom
    koes = pv2koe(truth_states_chief(:,i-1), MU);
    % ignoring J2 effects for now (likely ever)
    new_nu = keplerian_propagation( ...
        koes(6), ...
        e, ...
        a0, ...
        dt_now);
    [r_true, v_true] = koe2pv([koes(1:5).';new_nu], MU);
    truth_states_chief(:,i) = [r_true; v_true];
    
    %% Kalman Filter Time Update
    % state update
    [~, propagated_state] = ode45(@newton3d, [0, dt_now], estimated_states_chief(:, i-1), opts);
    state_minus = propagated_state(end, :).';

    pre_measurement_state_estimate(:, i) = state_minus;
    
    % covariance update
    r_vec = estimated_states_chief(1:3, i-1);
    jac = newtonian_jacobian(r_vec, MU);
    stm = stm_from_jacobian(jac, dt_now);
    
    % TODO: add process noise
    Pk_minus = stm * Ps(:,:,i-1) * stm.';

    %%  Kalman Filter Measurement Update
    % state update based off measurement
    measurement = measurement_model(truth_states_chief(:,i), measurement_noise);
    
    H = eye(6); % because dh/dx = I if we're only modelling noise?
    K = 0.01 * eye(6); % because I'm not sure what else to put here.
    
    % FIXME: measurement model feels sus
    state_plus = state_minus - K*(measurement_model(state_minus, measurement_noise) - measurement);
    Pk_plus = Pk_minus + K*H*Pk_minus*H.'*K.';
    
    % pack it all up
    estimated_states_chief(:,i) = state_plus;
    Ps(:,:,i) = Pk_plus;
end

%% Plot
figure

plot_t = ts / Tp;
subplot(4,2,1);
hold on
plot(plot_t, truth_states_chief(1, :), "DisplayName", "Truth")
plot(plot_t, estimated_states_chief(1, :), "DisplayName", "Estimated")
plot(plot_t, pre_measurement_state_estimate(1, :), "DisplayName", "Pre-Measurement Estimate")
ylabel("X Position [m]")
legend;
grid on;

subplot(4,2,3);
hold on
plot(plot_t, truth_states_chief(2, :), "DisplayName", "Truth")
plot(plot_t, estimated_states_chief(2, :), "DisplayName", "Estimated")
plot(plot_t, pre_measurement_state_estimate(2, :), "DisplayName", "Pre-Measurement Estimate")
ylabel("Y Position [m]")
legend;
grid on;

subplot(4,2,5);
hold on
plot(plot_t, truth_states_chief(3, :), "DisplayName", "Truth")
plot(plot_t, estimated_states_chief(3, :), "DisplayName", "Estimated")
plot(plot_t, pre_measurement_state_estimate(3, :), "DisplayName", "Pre-Measurement Estimate")

ylabel("Z Position [m]")
xlabel("Time [Orbital Periods]")
legend;
grid on;

subplot(4,2,7);
hold on
plot(plot_t, vecnorm(estimated_states_chief(1:3, :) - truth_states_chief(1:3, :)), "DisplayName", "Velocity Error")
ylabel("Position Error Magnitude [m]")
xlabel("Time [Orbital Periods]")
yscale log


subplot(4,2,2);
hold on
plot(plot_t, truth_states_chief(4, :), "DisplayName", "Truth")
plot(plot_t, estimated_states_chief(4, :), "DisplayName", "Estimated")
plot(plot_t, pre_measurement_state_estimate(4, :), "DisplayName", "Pre-Measurement Estimate")
ylabel("X Velocity [m/s]")
legend;
grid on;

subplot(4,2,4);
hold on
plot(plot_t, truth_states_chief(5, :), "DisplayName", "Truth")
plot(plot_t, estimated_states_chief(5, :), "DisplayName", "Estimated")
plot(plot_t, pre_measurement_state_estimate(5, :), "DisplayName", "Pre-Measurement Estimate")

ylabel("Y Velocity [m/s]")
legend;
grid on;

subplot(4,2,6);
hold on
plot(plot_t, truth_states_chief(6, :), "DisplayName", "Truth")
plot(plot_t, estimated_states_chief(6, :), "DisplayName", "Estimated")
plot(plot_t, pre_measurement_state_estimate(6, :), "DisplayName", "Pre-Measurement Estimate")
ylabel("Z Velocity [m/s]")
legend;
grid on;

subplot(4,2,8);
hold on
plot(plot_t, vecnorm(estimated_states_chief(4:6, :) - truth_states_chief(4:6, :)), "DisplayName", "Velocity Error")
ylabel("Velocity Error Magnitude [m/s]")
xlabel("Time [Orbital Periods]")
yscale log