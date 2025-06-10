close all
clear all
% Kalman Filter time, baby!
% constants
MU = 3.9860043550702260E+14; % m^3/s^2
J2 = 0.0010826358191967; % Earth's J2 coefficient
rE = 6.378136300e6; % Earth radius [m]

opts = odeset("RelTol",1e-12,"AbsTol",1e-16);


% initial absolute orbital elements
% initial Keplerian Orbital elements of Chief spacecraft (MMS-1):
ac = 9.758384355377886E+07; % semi-major axis, meters
ec = 0.1; % eccentricity
ic = 0.1; % inclination, degrees
raanc = 0; % Right Ascencion of the Ascending Node (aka longitude of ascending node aka omega), degrees
argpc = 0; % Argument of periapsis (aka Argument of Perifocus), degrees
nuc = 1; % True Anomaly, degrees

AOE_chief = [ac, ec, ic, raanc, argpc, nuc].';

[rc, vc] = koe2pv([ac, ec, ic, raanc, argpc, nuc], MU);

% initial normalized ROEs of the deputy spacecraft
ROEs = [0, 1000, 0, 0, 0, 0].' / ac; % unitless
AOE_deputy = roe2aoe(AOE_chief, ROEs);
[rd, vd] = koe2pv(AOE_deputy, MU);

initial_states_chief = [rc; vc];
initial_states_deputy = [rd; vd];


Tp = 2*pi / mean_motion(MU, ac);
num_steps = 1e3;
ts = linspace(100, 3*Tp, num_steps);

% setup state storage arrays
truth_states_chief = zeros(6, length(ts));
truth_states_deputy = zeros(6, length(ts));
measurements_chief = zeros(6, length(ts));
measurements_deputy = zeros(6, length(ts));
pre_measurement_state_estimate = zeros(6, length(ts));
truth_states_chief(:,1) = initial_states_chief;
truth_states_deputy(:,1) = initial_states_deputy;

estimated_states = zeros(6, length(ts));

% initial_offsets_chief = [10; -10; 10; 0.1; 0.1; 0.1];
initial_offsets_chief = [1; -1; 1; 0.01; 0.01; 0.01];
initial_offsets_chief = [0;0;0;0;0;0];
initial_offsets_deputy = -1 * initial_offsets_chief;

measurement_noise = [5; 5; 5; 5; 5; 5];

gps_measurement_chief = initial_states_chief + initial_offsets_chief;
gps_measurement_deputy = initial_states_deputy + initial_offsets_deputy;

initial_state_estimate = pv2roe(gps_measurement_chief, gps_measurement_deputy);

initial_state_estimate = ROEs + [100,100,100,100,100,100].' / ac;
estimated_states(:,1) = initial_state_estimate;
pre_measurement_state_estimate(:,1) = initial_state_estimate;
measurements(:,1) = initial_state_estimate;

initial_error = initial_state_estimate - ROEs;

Ps = zeros(6,6,length(ts));
P0 = diag(abs(initial_error));
Ps(:,:,1) = P0;

% from J2_STM_validation, the dynamics uncertainties are:
uncertainties = [5, 1, 2, 2, 0.5, 0.5].'; % meters
Q = uncertainties / ac;

% TODO(Office Hours): how go to ROE space properly? Answer: GOTO D'Amico's
% PhD, use linear mapping.
R = diag(measurement_noise.^3 / ac);

H = eye(6); % because dh/dx = I if we're only modelling noise?

for i=2:length(ts)
    %% Propagate Simulated Truth State forward in time
    dt_now = ts(i) - ts(i-1);
    [~, propagated_state_chief] = ode45(@newton3d_J2, [0, dt_now], truth_states_chief(:, i-1), opts);
    [~, propagated_state_deputy] = ode45(@newton3d_J2, [0, dt_now], truth_states_deputy(:, i-1), opts);

    truth_states_chief(:,i) = propagated_state_chief(end, :).';
    truth_states_deputy(:,i) = propagated_state_deputy(end, :).';

    %% Kalman Filter Time Update
    % state update
    % TODO: include control inputs
    previous_estimate = estimated_states(:, i-1);
    stm = chernick_J2_stm(pv2koe(truth_states_chief(:, i), MU), dt_now, rE, MU, J2);
    state_minus = stm * previous_estimate;
    pre_measurement_state_estimate(:, i) = state_minus;

    % covariance update
    % TODO: include control inputs
    Pk_minus = stm * Ps(:,:,i-1) * stm.' + Q;

    %%  Kalman Filter Measurement Update
    % state update based off measurement
    % [state_measurement, chief_measurement, deputy_measurement] = measurement_model(truth_states_chief(:,i), truth_states_deputy(:,i), measurement_noise);
    state_measurement = pv2roe(truth_states_chief(:,i), truth_states_deputy(:,i)) + measurement_noise .* randn(6,1) / ac;

    measurements(:,i) = state_measurement;

    K = Pk_minus * H.' / (H*Pk_minus * H.' + R);
    % modelled_measurement = state_minus + measurement_noise.* randn(6,1) .* [1, 1, 1, 0.1, 0.1, 0.1].';
    modelled_measurement = state_minus + measurement_noise.* randn(6,1) / ac;
    state_plus = state_minus + K*(state_measurement - modelled_measurement);

    % Joseph Formulation
    Pk_plus = (eye(6) - K*H) * Pk_minus * (eye(6) - K*H).' + K*R*K.';

    % pack it all up
    estimated_states(:,i) = state_plus;
    Ps(:,:,i) = Pk_plus;
end

%% Post-Processing
P_mags = zeros(6, length(ts));
true_roes = zeros(6, length(ts));
for k=1:length(ts)
    P_mags(:,k) = diag(Ps(:,:,k));
    true_roes(:,k) = pv2roe(truth_states_chief(:,i), truth_states_deputy(:,i));
end

%% Plot
figure
p=1;
plot_t = ts / Tp;
subplot(2,3,1);
hold on
plot(plot_t, ac * true_roes(p, :), "DisplayName", "Truth")
% plot(plot_t, ac * measurements(p, :), "DisplayName", "Measured")
plot(plot_t, ac * estimated_states(p, :), "DisplayName", "Estimated")
% plot(plot_t, ac * pre_measurement_state_estimate(p, :), "DisplayName", "Pre-Measurement Estimate")
ylabel("$a \delta a$ [m]", Interpreter='latex')
legend;
grid on;

p=2;
subplot(2,3, 2);
hold on
plot(plot_t, ac * true_roes(p, :), "DisplayName", "Truth")
% plot(plot_t, ac * measurements(p, :), "DisplayName", "Measured")
plot(plot_t, ac * estimated_states(2, :), "DisplayName", "Estimated")
% plot(plot_t, ac * pre_measurement_state_estimate(2, :), "DisplayName", "Pre-Measurement Estimate")
ylabel("$a \delta \lambda$ [m]", Interpreter='latex')
grid on;

p=3;
subplot(2,3,3);
hold on
plot(plot_t, ac * true_roes(p, :), "DisplayName", "Truth")
% plot(plot_t, ac * measurements(p, :), "DisplayName", "Measured")
plot(plot_t, ac * estimated_states(p, :), "DisplayName", "Estimated")
% plot(plot_t, ac * pre_measurement_state_estimate(p, :), "DisplayName", "Pre-Measurement Estimate")
ylabel("$a \delta e_x$ [m]", Interpreter='latex')
grid on;

p=4;
subplot(2,3,4);
hold on
plot(plot_t, ac * true_roes(p, :), "DisplayName", "Truth")
% plot(plot_t, ac * measurements(p, :), "DisplayName", "Measured")
plot(plot_t, ac * estimated_states(p, :), "DisplayName", "Estimated")
% plot(plot_t, ac * pre_measurement_state_estimate(p, :), "DisplayName", "Pre-Measurement Estimate")

ylabel("$a \delta e_y$ [m]", Interpreter='latex')


p=5;
subplot(2,3,5);
hold on
plot(plot_t, ac * true_roes(p, :), "DisplayName", "Truth")
% plot(plot_t, ac * measurements(p, :), "DisplayName", "Measured")
plot(plot_t, ac * estimated_states(p, :), "DisplayName", "Estimated")
% plot(plot_t, ac * pre_measurement_state_estimate(p, :), "DisplayName", "Pre-Measurement Estimate")

ylabel("$a \delta i_x$ [m]", Interpreter='latex')


p=6;
subplot(2,3,6);
hold on
plot(plot_t, ac * true_roes(p, :), "DisplayName", "Truth")
% plot(plot_t, ac * measurements(p, :), "DisplayName", "Measured")
plot(plot_t, ac * estimated_states(p, :), "DisplayName", "Estimated")
% plot(plot_t, ac * pre_measurement_state_estimate(p, :), "DisplayName", "Pre-Measurement Estimate")

ylabel("$a \delta i_y$ [m]", Interpreter='latex')
xlabel("Time [Orbital Periods]")

figure

p=1;
plot_t = ts / Tp;
subplot(2,3,1);
hold on
% plot(plot_t, (ac * (measurements(p, :) - true_roes(p, :))), "DisplayName", "Measurement Error")
plot(plot_t, (ac * (estimated_states(p, :) - true_roes(p, :))), "DisplayName", "Estimate Error")
% plot(plot_t, abs(ac * (pre_measurement_state_estimate(p, :)- true_roes(p, :))), "DisplayName", "Pre-Measurement Estimate Error")
fill([plot_t, fliplr(plot_t)], ac * [P_mags(p,:), fliplr(-P_mags(p,:))], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName','Covariance Magnitude')
ylabel("$\Delta a \delta a$ [m]", Interpreter='latex')

legend;
grid on;

p=2;
subplot(2,3,2);
hold on
% plot(plot_t, (ac * (measurements(p, :) - true_roes(p, :))), "DisplayName", "Measurement Error")
plot(plot_t, (ac * (estimated_states(p, :) - true_roes(p, :))), "DisplayName", "Estimate Error")
% plot(plot_t, abs(ac * (pre_measurement_state_estimate(p, :)- true_roes(p, :))), "DisplayName", "Pre-Measurement Estimate Error")
fill([plot_t, fliplr(plot_t)], ac * [P_mags(p,:), fliplr(-P_mags(p,:))], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
ylabel("$\Delta a \delta \lambda$ [m]", Interpreter='latex')

grid on;

p=3;
subplot(2,3,3);
hold on
% plot(plot_t, (ac * (measurements(p, :) - true_roes(p, :))), "DisplayName", "Measurement Error")
plot(plot_t, (ac * (estimated_states(p, :) - true_roes(p, :))), "DisplayName", "Estimate Error")
% plot(plot_t, abs(ac * (pre_measurement_state_estimate(p, :)- true_roes(p, :))), "DisplayName", "Pre-Measurement Estimate Error")
fill([plot_t, fliplr(plot_t)], ac * [P_mags(p,:), fliplr(-P_mags(p,:))], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
ylabel("$\Delta a \delta e_x$ [m]", Interpreter='latex')

grid on;

p=4;
subplot(2,3,4);
hold on
% plot(plot_t, (ac * (measurements(p, :) - true_roes(p, :))), "DisplayName", "Measurement Error")
plot(plot_t, (ac * (estimated_states(p, :) - true_roes(p, :))), "DisplayName", "Estimate Error")
% plot(plot_t, abs(ac * (pre_measurement_state_estimate(p, :)- true_roes(p, :))), "DisplayName", "Pre-Measurement Estimate Error")
fill([plot_t, fliplr(plot_t)], ac * [P_mags(p,:), fliplr(-P_mags(p,:))], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
ylabel("$\Delta a \delta e_y$ [m]", Interpreter='latex')



p=5;
subplot(2,3,5);
hold on
% plot(plot_t, (ac * (measurements(p, :) - true_roes(p, :))), "DisplayName", "Measurement Error")
plot(plot_t, (ac * (estimated_states(p, :) - true_roes(p, :))), "DisplayName", "Estimate Error")
% plot(plot_t, abs(ac * (pre_measurement_state_estimate(p, :)- true_roes(p, :))), "DisplayName", "Pre-Measurement Estimate Error")
fill([plot_t, fliplr(plot_t)], ac * [P_mags(p,:), fliplr(-P_mags(p,:))], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
ylabel("$\Delta a \delta i_x$ [m]", Interpreter='latex')



p=6;
subplot(2,3,6);
hold on
% plot(plot_t, (ac * (measurements(p, :) - true_roes(p, :))), "DisplayName", "Measurement Error")
plot(plot_t, (ac * (estimated_states(p, :) - true_roes(p, :))), "DisplayName", "Estimate Error")
% plot(plot_t, abs(ac * (pre_measurement_state_estimate(p, :)- true_roes(p, :))), "DisplayName", "Pre-Measurement Estimate Error")
fill([plot_t, fliplr(plot_t)], ac * [P_mags(p,:), fliplr(-P_mags(p,:))], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')

ylabel("$\Delta a \delta i_y$ [m]", Interpreter='latex')
xscale log

xlabel("Time [Orbital Periods]")