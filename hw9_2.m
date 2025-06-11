% HW9
% Oh yeah, it's all coming together
close all; clear; clc;

% Constants
J2 = 0.0010826358191967; % Earth's J2 coefficient
mu = 3.986004415e14; % Earth gravitational parameter [m^3/s^2]
rE = 6.378136300e6; % Earth radius [m]

% Plotting
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesFontSize',18);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


% Integration
opts = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances

% %%%%%%%%%%%
%% Problem 2c
% %%%%%%%%%%%

% Kalman Filter time, baby!

% initial absolute orbital elements
% initial Keplerian Orbital elements of Chief spacecraft (MMS-1):
ac = 9.758384355377886E+07; % semi-major axis, meters
ec = 0.001; % eccentricity
ic = 0.1; % inclination, degrees
raanc = 0; % Right Ascencion of the Ascending Node (aka longitude of ascending node aka omega), degrees
argpc = 0; % Argument of periapsis (aka Argument of Perifocus), degrees
nuc = 1; % True Anomaly, degrees

n = sqrt(mu/ac^3); % mean motion


AOE_chief = [ac, ec, ic, raanc, argpc, nuc].';

[rc, vc] = koe2pv([ac, ec, ic, raanc, argpc, nuc], mu);

% initial normalized ROEs of the deputy spacecraft
ROEs = [0, 1000, 0, 0, 0, 0].' / ac; % unitless
AOE_deputy = roe2aoe(AOE_chief, ROEs);
[rd, vd] = koe2pv(AOE_deputy, mu);

initial_states_chief = [rc; vc];
initial_states_deputy = [rd; vd];


Tp = 2*pi / mean_motion(mu, ac);
num_steps = 1e3;
ts = linspace(100, 3*Tp, num_steps);

% setup state storage arrays
truth_states_chief = [];
truth_states_deputy = [];
pre_measurement_state_estimate = [];
truth_states_chief(:,1) = initial_states_chief;
truth_states_deputy(:,1) = initial_states_deputy;

estimated_states = [];

% initial_offsets_chief = [10; -10; 10; 0.1; 0.1; 0.1];
initial_offsets_chief = [1; -1; 1; 0.01; 0.01; 0.01];
initial_offsets_chief = [0;0;0;0;0;0];
initial_offsets_deputy = -1 * initial_offsets_chief;

measurement_noise = [5; 5; 5; 5; 5; 5];

gps_measurement_chief = initial_states_chief + initial_offsets_chief;
gps_measurement_deputy = initial_states_deputy + initial_offsets_deputy;

% initial_state_estimate = pv2roe(gps_measurement_chief, gps_measurement_deputy);

initial_state_estimate = ROEs + [100,100,100,100,100,100].' / ac;
estimated_states(:,1) = initial_state_estimate;
pre_measurement_state_estimate(:,1) = initial_state_estimate;
measurements(:,1) = initial_state_estimate;

initial_error = initial_state_estimate - ROEs;

Ps = zeros(6,6,1);
P0 = diag(abs(initial_error));
Ps(:,:,1) = P0;

% from J2_STM_validation, the dynamics uncertainties are:
uncertainties = [5, 1, 2, 2, 0.5, 0.5].'; % meters
Q = uncertainties / ac;

% TODO(Office Hours): how go to ROE space properly? Answer: GOTO D'Amico's
% PhD, use linear mapping.
R = diag(measurement_noise.^3 / ac);

H = eye(6); % because dh/dx = I if we're only modelling noise

DESIRED_ROES = [1000, 1000, 0, 0, 1000, 0].';
DT = 5*Tp;

dvs = [];
% let's say we have 5 time regions of stationkeeping
for k=1:5
    % allow a few orbits for the EKF to converge, then start stationkeeping
    if mod(k, 2) == 0
        % plan maneuvers based off of current state estimate
        % Define the functions you'd like to use for control input matrix and STM
        stm = @(chief_oe, t) chernick_J2_stm(chief_oe, t, rE, mu, J2);
        control_input_matrix = @(chief_oe) chernick_control_matrix(chief_oe, mu);
        koe = pv2koe(propagated_state_chief(end, :).', mu);
        [t_maneuvers, maneuvers, total_cost] = impulsive_control(koe, koe(1)*state_plus, DESIRED_ROES, DT, stm, control_input_matrix, rE, mu, J2);
    else
        t_maneuvers = [1];
        maneuvers = [0; 0; 0;];
    end
    propagation_times_start = [0, t_maneuvers];
    propagation_times_end = [t_maneuvers, DT];
    for j = 1:length(t_maneuvers)+1
        ts = linspace(propagation_times_start(j), propagation_times_end(j));
        for i=2:length(ts)
            %% Propagate Simulated Truth State forward in time
            dt_now = ts(i) - ts(i-1);
            [~, propagated_state_chief] = ode45(@newton3d_J2, [0, dt_now], truth_states_chief(:, end), opts);
            [~, propagated_state_deputy] = ode45(@newton3d_J2, [0, dt_now], truth_states_deputy(:, end), opts);

            truth_states_chief = [truth_states_chief propagated_state_chief(end, :).'];
            truth_states_deputy = [truth_states_deputy propagated_state_deputy(end, :).'];

            %% Kalman Filter Time Update
            % state update
            previous_estimate = estimated_states(:, end);
            phi = chernick_J2_stm(pv2koe(truth_states_chief(:, end), mu), dt_now, rE, mu, J2);
            state_minus = phi * previous_estimate;
            pre_measurement_state_estimate =[pre_measurement_state_estimate, state_minus];

            % covariance update
            Pk_minus = phi * Ps(:,:,end) * phi.' + Q;

            %%  Kalman Filter Measurement Update
            % state update based off measurement
            % [state_measurement, chief_measurement, deputy_measurement] = measurement_model(truth_states_chief(:,i), truth_states_deputy(:,i), measurement_noise);
            state_measurement = pv2roe(truth_states_chief(:,end), truth_states_deputy(:,end)) + measurement_noise .* randn(6,1) / ac;

            measurements = [measurements, state_measurement];

            K = Pk_minus * H.' / (H*Pk_minus * H.' + R);
            % modelled_measurement = state_minus + measurement_noise.* randn(6,1) .* [1, 1, 1, 0.1, 0.1, 0.1].';
            modelled_measurement = state_minus + measurement_noise.* randn(6,1) / ac;
            state_plus = state_minus + K*(state_measurement - modelled_measurement);

            % Joseph Formulation
            Pk_plus = (eye(6) - K*H) * Pk_minus * (eye(6) - K*H).' + K*R*K.';

            % pack it all up
            estimated_states = [estimated_states, state_plus];
            Ps = cat(3, Ps, Pk_plus);
        end

        % apply maneuver dv
        if j <= length(t_maneuvers)
            last_state = truth_states_deputy(:, end);
            r_last = last_state(1:3);
            v_last = last_state(4:6);
            dv_RTN = maneuvers(:, j);
            disp(dv_RTN);
            % rotate vector:
            dv_ECI = 1e0*rtn2eci(r_last, v_last, dv_RTN);
            state_override = [0; 0; 0; dv_ECI];
            
            truth_states_deputy(:, end) = truth_states_deputy(:, end) + state_override;
            disp(truth_states_deputy(:, end)-last_state);
        else
            % truth_states_deputy(:, end) = [];
            % truth_states_chief(:, end) = [];
        end
        %%
        % TODO: Kalman filter update with control input!
        

        % dv accounting:
        if j > 1
            prior_dv = dvs(end);
            dv_mag = norm(dv_RTN);
        else
            dv_mag = 0;
            prior_dv = 0;
        end

        dvs = [dvs, (prior_dv+dv_mag)*ones(1,length(ts))];
    end
end

%% Post-Processing
P_mags = zeros(6, length(Ps));
true_roes = zeros(6, length(truth_states_chief));
for k=1:length(Ps)
    P_mags(:,k) = diag(Ps(:,:,k));
    true_roes(:,k) = pv2roe(truth_states_chief(:,k), truth_states_deputy(:,k));
end
true_roes(6,1) = 0;

%% Plot
figure
p=1;
% plot_t = ts / Tp;
plot_t = 1:length(truth_states_chief);
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
% plot_t = ts / Tp;
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