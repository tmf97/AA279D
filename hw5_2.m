close all
clear truth_ephem
clear prop_ephem
% validate propagator:
MU = 3.9860043550702260E+14; % m^3/s^2
opts = odeset("RelTol",1e-13,"AbsTol",1e-15);

% for purposes of this analysis, have the deputy be at +10km in T at apogee
% read truth ephemerides:
% one column vector for all state.
[mms1_ephem, ts] = read_horizons('mms-1.txt');
[mms2_ephem, ~]  = read_horizons('mms-2.txt');

% load initial conditions
aoe_chief_i = mms1_ephem(:,1);
aoe_deputy_i = mms2_ephem(:,1);


% convert osculating to mean
maoe_chief_i = osc2mean(aoe_chief_i, 1);
maoe_deputy_i = osc2mean(aoe_deputy_i, 1);
a0 = maoe_chief_i(1);
% propagate the true anomaly of the chief out to apogee
nu_apogee = pi;
maoe_chief_f = [maoe_chief_i(1:end-1); nu_apogee];

% convert keplerian to cartesian
[pos_chief_f_eci, vel_chief_f_eci]  = koe2pv(maoe_chief_f, MU);

% convert position in RTN frame to cartesian
deputy_pos_RTN = [0; 0; 10000]; % meters
deputy_pos_ECI = pos_chief_f_eci + rtn2eci(pos_chief_f_eci, vel_chief_f_eci, deputy_pos_RTN);

% convert desired orbital elements
% assume deputy velocity = chief velocity
aoe_deputy_desired = pv2koe([deputy_pos_ECI; vel_chief_f_eci], MU);


% delta alpha bar final
mean_roes_desired = quasi_nonsingular_roe(maoe_chief_f, aoe_deputy_desired);

% delta alpha bar initial
mean_roes_initial = quasi_nonsingular_roe(maoe_chief_i, maoe_deputy_i);

% find time at apogee relative to initial
ti = 0;
tf = orbit_timing(maoe_chief_i, nu_apogee, MU);
Tp = 2*pi/mean_motion(MU, a0);
tf = tf + Tp; % workaround


% compute STM 
stm_fi = STM_J2(maoe_chief_i, ti, tf, MU);

% compute pseudostate
w = mean_roes_desired - stm_fi * mean_roes_initial;

% define manuevers a priori
maneuver_nus = [deg2rad(150), deg2rad(220), deg2rad(0), deg2rad(45)]; % keep manuevers at lower altitude to not interfere with science, if possible
maneuver_ts = zeros(size(maneuver_nus));
for i=1:length(maneuver_nus)
    nu = maneuver_nus(i);
    t = orbit_timing(maoe_chief_i, nu, MU);
    maneuver_ts(i) = t;
end
% compute M matricies for each maneuver and concatenate

Ms = [];

for k=1:length(maneuver_nus)
    nu = maneuver_nus(k);
    t = maneuver_ts(k);
    stm_fk = STM_kep(maoe_chief_i, t, tf, MU);
    Gamma_k = control_input_matrix(maoe_chief_i, nu, MU);
    M = stm_fk * Gamma_k;
    Ms = [Ms, M];
end

% compute maneuvers from pseudoinverse of M
dvs = lsqminnorm(Ms,w);

% each column is the dv of the maneuver in RTN
maneuver_dvs = reshape(dvs, [3, length(maneuver_ts)]);

maneuver_names = string(rad2deg(maneuver_nus)) + " deg";
disp("Least-Squares Manevers")
T = array2table(maneuver_dvs, 'VariableNames', maneuver_names, 'RowName', {'dv_R', 'dv_T', 'dv_N'});
disp(T);


% propagate orbit to each maneuver, apply dv manually
[rc0, vc0] = koe2pv(aoe_chief_i, MU);
[rd0, vd0] = koe2pv(aoe_deputy_i, MU);
states_chief = [rc0; vc0];
states_deputy = [rd0; vd0];
propagation_times_start = [0, maneuver_ts];
propagation_times_end = [maneuver_ts, tf];
times = [];
for i = 1:length(maneuver_ts)+1
    t_start = propagation_times_start(i);
    t_end = propagation_times_end(i);
    ts = linspace(t_start, t_end, 1e4);
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
    if i <= length(maneuver_ts)
        last_state = states_deputy(:, end);
        r_last = last_state(1:3);
        v_last = last_state(4:6);
        dv_RTN = maneuver_dvs(:, i);
        % rotate vector:
        rhat = r_last/norm(r_last);
        h0 = cross(r_last, v_last);
        nhat = h0/norm(h0);
        that = cross(nhat, rhat);
        C_RTN_to_ECI = [rhat, that, nhat];
        
        dv_ECI = C_RTN_to_ECI*dv_RTN;
        state_override = [0;0;0; dv_ECI];
        states_deputy(:, end) = states_deputy(:, end) + state_override;
    else
        states_deputy(:, end) = [];
        states_chief(:, end) = [];
    end
end
% plot everything in RTN of the chief:
mms2_rtns = zeros(3, length(states_chief));
mms2_rtn_vs = zeros(3, length(states_chief));
mean_roes = zeros(6, length(states_chief));
for i=1:length(states_chief)
    r_chief = states_chief(1:3,i);
    v_chief = states_chief(4:6,i);

    r_deputy = states_deputy(1:3,i);
    v_deputy = states_deputy(4:6,i);
    [mms2_rtns(:,i), mms2_rtn_vs(:,i)] = pvs2rtn(r_chief, v_chief, r_deputy, v_deputy);
    % compute KOEs for both chief and deputy
    koe_chief = pv2koe([r_chief; v_chief], MU);
    koe_deputy = pv2koe([r_deputy; v_deputy], MU);

    moe_chief = osc2mean(koe_chief, 1);
    moe_deputy = osc2mean(koe_deputy, 1);
    mean_roes(:, i) = quasi_nonsingular_roe(moe_chief, moe_deputy);
end

Rs2 = mms2_rtns(1, :);
Ts2 = mms2_rtns(2, :);
Ns2 = mms2_rtns(3, :);
times = times.';
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
axis square;



% 3D!
figure
grid on
hold on
axis equal
Re = 6378000; % m
[X,Y,Z] = sphere;
surf(X*Re, Y*Re, Z*Re, 'DisplayName', 'Earth')
plot3(states_chief(1,:), states_chief(2,:), states_chief(3,:))
plot3(states_deputy(1,:), states_deputy(2,:), states_deputy(3,:))


figure 
hold on
plot(times, vecnorm(states_chief(1:3,:)));
plot(times, vecnorm(states_deputy(1:3,:)));
ylabel("[m]")


% plot ROEs
figure
subplot(1,3,1)
hold on
plot(a0*mean_roes(2, :), a0*mean_roes(1, :), 'DisplayName', 'Trajectory')
scatter(mean_roes_desired(2), mean_roes_desired(1), 'DisplayName', 'Target')
scatter(a0*mean_roes(2, 1), a0*mean_roes(1, 1), 'DisplayName','Start');
xlabel("$a \delta \lambda$ (m)")
ylabel("$a \delta a$ (m)")
grid on;
axis square;

hold off
legend


subplot(1,3,2)
hold on
plot(a0*mean_roes(3, :), a0*mean_roes(4, :))
scatter(mean_roes_desired(3), mean_roes_desired(4))
scatter(a0*mean_roes(3, 1), a0*mean_roes(4, 1));
xlabel("$a \delta e_x$ (m)")
ylabel("$a \delta e_y$ (m)")
grid on;
axis square;

hold off

subplot(1,3,3)
hold on
plot(a0*mean_roes(5, :), a0*mean_roes(6, :))
scatter(mean_roes_desired(5), mean_roes_desired(6))
scatter(a0*mean_roes(5, 1), a0*mean_roes(6, 1));
xlabel("$a \delta i_x$ (m)")
ylabel("$a \delta i_y$ (m)")
grid on;
axis square;

hold off

