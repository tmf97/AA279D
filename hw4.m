close all
% chief AOEs
AOE_chief = [63781000, 0.001, deg2rad(30.11885), deg2rad(356.7045), deg2rad(158.2343), deg2rad(0.1)];
ROEs = [0, 100, 50, 100, 30, 200]; % meters

AOE_deputy = roe2aoe(AOE_chief, ROEs);

MU = 3.9860043550702260E+14; % m^3/s^2

disp_aoe_chief = [
    AOE_chief(1)/1000;
    AOE_chief(2);
    rad2deg(AOE_chief(3));
    rad2deg(AOE_chief(4));
    rad2deg(AOE_chief(5));
    rad2deg(AOE_chief(6));
    ];

disp_aoe_deputy = [
    AOE_deputy(1)/1000;
    AOE_deputy(2);
    rad2deg(AOE_deputy(3));
    rad2deg(AOE_deputy(4));
    rad2deg(AOE_deputy(5));
    rad2deg(AOE_deputy(6));
    ];

disp('Initial Orbital Elements:')
T = array2table([disp_aoe_chief.'; disp_aoe_deputy.'; (disp_aoe_chief-disp_aoe_deputy).'], 'VariableNames', {'a (km)', 'e', 'i (deg)', 'Omega (deg)', 'w (deg)', 'f (deg)'}, 'RowName', {'Chief','Deputy', 'Difference'});
disp(T);



[rc, vc] = koe2pv(AOE_chief, MU);
[rd, vd] = koe2pv(AOE_deputy, MU);
sc = [rc;vc];
sd = [rd;vd];

disp('Initial Positions and Velocities:')
T = array2table([sc.'; sd.'; (sd-sc).'], 'VariableNames', {'r_x (m)', 'r_y (m)', 'r_z (m)', 'v_x (m/s)', 'v_y (m/s)', 'v_z (m/s)'}, 'RowName', {'Chief','Deputy', 'Difference'});
disp(T);

ac = AOE_chief(1);
Tp = 2*pi*sqrt(ac^3 / MU);

ts = linspace(0, 15*Tp, 1e4);

opts = odeset("RelTol",1e-13,"AbsTol",1e-16);



[~, states_chief] = ode45(@newton3d, ts, sc, opts);
[~, states_deputy] = ode45(@newton3d, ts, sd, opts);
kep_states_chief = zeros(size(states_chief));
kep_states_deputy = zeros(size(states_deputy));
roes = zeros(size(states_deputy));
mean_oes_chief = zeros(size(states_chief));
mean_oes_deputy = zeros(size(states_deputy));
mean_roes = zeros(size(states_deputy));
rtns = zeros(size(states_deputy));

for i=1:length(ts)
    state_chief = states_chief(i,:);
    state_deputy = states_deputy(i,:);
    kep_chief = pv2koe(state_chief, MU);
    kep_deputy = pv2koe(state_deputy, MU);
    kep_states_chief(i,:) = kep_chief;
    kep_states_deputy(i,:) = kep_deputy;
    roe = quasi_nonsingular_roe(kep_chief, kep_deputy);
    roes(i,:) = roe;

    % mean orbital elements
    mean_oe_chief = osc2mean(kep_chief, false);
    mean_oe_deputy = osc2mean(kep_deputy, false); 
    
    mean_roe = quasi_nonsingular_roe(mean_oe_chief, mean_oe_deputy);

    mean_oes_chief(i,:) = mean_oe_chief;
    mean_oes_deputy(i,:) = mean_oe_deputy;
    mean_roes(i,:) = mean_roe;

    % RTN stuff
    r0 = state_chief(1:3);
    v0 = state_chief(4:6);
    r1 = state_deputy(1:3);
    v1 = state_deputy(4:6);
    ac = kep_chief(1);
    [rtn, rtn_dot] = pvs2rtn(r0, v0,r1, v1, MU, ac);
    rtns(i,:) = [rtn; rtn_dot];
end

[~, states_chief_j2] = ode45(@newton3d_J2, ts, sc, opts);
[~, states_deputy_j2] = ode45(@newton3d_J2, ts, sd, opts);
kep_states_chief_j2 = zeros(size(states_chief_j2));
kep_states_deputy_j2 = zeros(size(states_deputy_j2));
roes_j2 = zeros(size(states_deputy_j2));
mean_oes_chief_j2 = zeros(size(states_chief_j2));
mean_oes_deputy_j2 = zeros(size(states_deputy_j2));
mean_roes_j2 = zeros(size(states_deputy_j2));
rtns_j2 = zeros(size(states_deputy_j2));

for i=1:length(ts)
    state_chief_j2= states_chief_j2(i,:);
    state_deputy_j2 = states_deputy_j2(i,:);
    kep_chief_j2 = pv2koe(state_chief_j2, MU);
    kep_deputy_j2 = pv2koe(state_deputy_j2, MU);
    kep_states_chief_j2(i,:) = kep_chief_j2;
    kep_states_deputy_j2(i,:) = kep_deputy_j2;
    roe_j2 = quasi_nonsingular_roe(kep_chief_j2, kep_deputy_j2);
    roes_j2(i,:) = roe_j2;

    % mean orbital elements
    mean_oe_chief_j2 = osc2mean(kep_chief_j2, true);
    mean_oe_deputy_j2 = osc2mean(kep_deputy_j2, true); 
    
    mean_roe_j2 = quasi_nonsingular_roe(mean_oe_chief_j2, mean_oe_deputy_j2);

    mean_oes_chief_j2(i,:) = mean_oe_chief_j2;
    mean_oes_deputy_j2(i,:) = mean_oe_deputy_j2;
    mean_roes_j2(i,:) = mean_roe_j2;

    % RTN stuff
    r0 = state_chief_j2(1:3);
    v0 = state_chief_j2(4:6);
    r1 = state_deputy_j2(1:3);
    v1 = state_deputy_j2(4:6);
    ac = kep_chief_j2(1);
    [rtn, rtn_dot] = pvs2rtn(r0, v0,r1, v1, MU, ac);
    rtns_j2(i,:) = [rtn; rtn_dot];
end


as = kep_states_chief(:,1);
es = kep_states_chief(:,2);
is = kep_states_chief(:,3);
Omegas = kep_states_chief(:,4);
ws = kep_states_chief(:,5);
nus = kep_states_chief(:,6);

mean_as = mean_oes_chief(:,1);
mean_es = mean_oes_chief(:,2);
mean_is = mean_oes_chief(:,3);
mean_Omegas = mean_oes_chief(:,4);
mean_ws = mean_oes_chief(:,5);
mean_nus = mean_oes_chief(:,6);

das = roes(:,1);
dlambdas = roes(:,2);
dexs = roes(:,3);
deys = roes(:,4);
dixs = roes(:,5);
diys = roes(:,6);

mean_das = mean_roes(:,1);
mean_dlambdas = mean_roes(:,2);
mean_dexs = mean_roes(:,3);
mean_deys = mean_roes(:,4);
mean_dixs = mean_roes(:,5);
mean_diys = mean_roes(:,6);

times = ts / Tp;

figure
subplot(6, 1, 1)
hold on
plot(times, as, 'DisplayName', 'Osculating');
plot(times, mean_as, 'DisplayName', 'Mean');
ylabel('Semi-major Axis (m)');
legend

subplot(6, 1, 2)
hold on
plot(times, es, 'DisplayName', 'Osculating');
plot(times, mean_es, 'DisplayName', 'Mean');
ylabel('Eccentricity')

subplot(6, 1, 3)
hold on
plot(times, is, 'DisplayName', 'Osculating');
plot(times, mean_is, 'DisplayName', 'Mean');
ylabel('Inclination (rad)')

subplot(6, 1, 4)
hold on
plot(times, Omegas, 'DisplayName', 'Osculating');
plot(times, mean_Omegas, 'DisplayName', 'Mean');
ylabel("RAAN (rad)")

subplot(6, 1, 5)
hold on
plot(times, ws, 'DisplayName', 'Osculating');
plot(times, mean_ws, 'DisplayName', 'Mean');
ylabel('Argument of Periapsis (rad)')
xlabel('Orbit number')


subplot(6, 1, 6)
hold on
plot(times, ws + nus - ws(1), 'DisplayName', 'Osculating');
plot(times, mean_ws + mean_nus - mean_ws(1), 'DisplayName', 'Mean');
ylabel('Arg of Latitude (rad)')
xlabel('Orbit number')


figure
subplot(6, 1, 1)
hold on
plot(times, das, 'DisplayName', 'Osculating');
plot(times, mean_das, 'DisplayName', 'Mean');
ylabel('$\delta a$', Interpreter='latex')
legend

subplot(6, 1, 2)
hold on
plot(times, dlambdas, 'DisplayName', 'Osculating');
plot(times, mean_dlambdas, 'DisplayName', 'Mean');
ylabel('$\delta \lambda$', Interpreter='latex')

subplot(6, 1, 3)
hold on
plot(times, dexs, 'DisplayName', 'Osculating');
plot(times, mean_dexs, 'DisplayName', 'Mean');
ylabel('$\delta e_x$', Interpreter='latex')

subplot(6, 1, 4)
hold on
plot(times, deys, 'DisplayName', 'Osculating');
plot(times, mean_deys, 'DisplayName', 'Mean');
ylabel('$\delta e_y$', Interpreter='latex')

subplot(6, 1, 5)
hold on
plot(times, dixs, 'DisplayName', 'Osculating');
plot(times, mean_dixs, 'DisplayName', 'Mean');
ylabel('$\delta i_x$', Interpreter='latex')

subplot(6, 1, 6)
hold on
plot(times, diys, 'DisplayName', 'Osculating');
plot(times, mean_diys, 'DisplayName', 'Mean');
ylabel('$\delta i_y$', Interpreter='latex')
xlabel('Orbit Number')


fontsize(12, 'points')


Rs = rtns(:, 1);
Ts = rtns(:, 2);
Ns = rtns(:, 3);

Rs_j2 = rtns_j2(:, 1);
Ts_j2 = rtns_j2(:, 2);
Ns_j2 = rtns_j2(:, 3);

figure

subplot(3,2,1);
hold on
plot(times, Rs);
plot(times, Rs_j2);
ylabel("R (m)")
grid on;


subplot(3,2,3);
hold on
plot(times, Ts);
plot(times, Ts_j2);
ylabel("T (m)")
grid on;


subplot(3,2,5);
hold on
plot(times, Ns);
plot(times, Ns_j2);
ylabel("N (m)")
xlabel("Orbital Periods")
grid on;


subplot(2,4,3)
hold on
plot(Ts, Rs)
plot(Ts_j2, Rs_j2)
xlabel("T (m)")
ylabel("R (m)")
grid on;
axis square;

subplot(2,4,4)
hold on
plot(Ns, Rs)
plot(Ns_j2, Rs_j2)
xlabel("N (m)")
ylabel("R (m)")
grid on;
axis square;

subplot(2,4,7)
hold on
plot(Ts, Ns)
plot(Ts_j2, Ns_j2)
xlabel("T (m)")
ylabel("N (m)")
grid on;
axis square;

subplot(2,4,8)
hold on
plot3(Rs, Ts, Ns)
plot3(Rs_j2, Ts_j2, Ns_j2)
xlabel("R (m)")
ylabel("T (m)")
zlabel("N (m)")
grid on;
axis square;


% plot difference between J2 and not in RTN

figure
subplot(3,2,1);
hold on
plot(times, Rs_j2 - Rs);
ylabel("\Delta R (m)")
grid on;


subplot(3,2,3);
hold on
plot(times, Ts_j2 - Ts);
ylabel("\Delta T (m)")
grid on;


subplot(3,2,5);
hold on
plot(times, Ns_j2 - Ns);
ylabel("\Delta N (m)")
xlabel("Orbital Periods")
grid on;


subplot(2,4,3)
hold on
plot(Ts_j2 - Ts, Rs_j2 - Rs)
xlabel("\Delta T (m)")
ylabel("\Delta R (m)")
grid on;
axis square;

subplot(2,4,4)
hold on
plot(Ns_j2 - Ns, Rs_j2 - Rs)
xlabel("\Delta N (m)")
ylabel("\Delta R (m)")
grid on;
axis square;

subplot(2,4,7)
hold on
plot(Ts_j2 - Ts, Ns_j2 - Ns)
xlabel("\Delta T (m)")
ylabel("\Delta N (m)")
grid on;
axis square;

subplot(2,4,8)
hold on
plot3(Rs_j2 - Rs, Ts_j2 - Ts, Ns_j2 - Ns)
xlabel("\Delta R (m)")
ylabel("\Delta T (m)")
zlabel("\Delta N (m)")
grid on;
axis square;


%%%%

das = ac * roes_j2(:,1);
dls = ac *roes_j2(:,2);
dexs = ac *roes_j2(:,3);
deys = ac *roes_j2(:,4);
dixs = ac *roes_j2(:,5);
diys = ac *roes_j2(:,6);

mdas = ac * mean_roes_j2(:,1);
mdls = ac *mean_roes_j2(:,2);
mdexs = ac *mean_roes_j2(:,3);
mdeys = ac *mean_roes_j2(:,4);
mdixs = ac *mean_roes_j2(:,5);
mdiys = ac *mean_roes_j2(:,6);


figure
subplot(1, 3, 1)
hold on
plot(dexs, deys)
plot(mdexs, mdeys)
xlabel('$a \delta e_x$ (m)', Interpreter='latex')
ylabel('$a \delta e_y$ (m)', Interpreter='latex')
axis square


subplot(1, 3, 2)
hold on
plot(dixs, diys, 'DisplayName','Osculating')
plot(mdixs, mdiys, 'DisplayName','Mean')
xlabel('$a \delta i_x$ (m)', Interpreter='latex')
ylabel('$a \delta i_y$ (m)', Interpreter='latex')
axis square
legend

subplot(1, 3, 3)
hold on
plot(dls, das)
plot(mdls, mdas)
xlabel('$a \delta \lambda$ (m)', Interpreter='latex')
ylabel('$a \delta a$ (m)', Interpreter='latex')

axis square
fontsize(11, 'points')



%%%%%%%%%%%%%%%%%%%%%%%%%

% set dix to 0 to remove drift in diy and dlambda

%%%%%%%%%%%%%%%%%%%%%%%%%%
ROEs = [0, 100, 50, 100, 0, 200]; % meters
AOE_deputy = roe2aoe(AOE_chief, ROEs);
[rd, vd] = koe2pv(AOE_deputy, MU);
sd = [rd;vd];
[~, states_deputy_J2] = ode45(@newton3d_J2, ts, sd, opts);
kep_states_deputy_j2 = zeros(size(states_deputy));
roes_j2 = zeros(size(states_deputy));
mean_oes_deputy_j2 = zeros(size(states_deputy));
mean_roes_j2 = zeros(size(states_deputy));


for i=1:length(ts)
    state_deputy = states_deputy_J2(i,:);
    kep_chief = kep_states_chief_j2(i,:);
    kep_deputy = pv2koe(state_deputy, MU);
    roe = quasi_nonsingular_roe(kep_chief, kep_deputy);
    roes_j2(i,:) = roe;

    % mean orbital elements
    mean_oe_chief = mean_oes_chief_j2(i,:);
    mean_oe_deputy = osc2mean(kep_deputy, true); 
    
    mean_roe = quasi_nonsingular_roe(mean_oe_chief, mean_oe_deputy);
    mean_oes_deputy_j2(i,:) = mean_oe_deputy;
    mean_roes_j2(i,:) = mean_roe;
end

das = ac * roes_j2(:,1);
dls = ac *roes_j2(:,2);
dexs = ac *roes_j2(:,3);
deys = ac *roes_j2(:,4);
dixs = ac *roes_j2(:,5);
diys = ac *roes_j2(:,6);

mdas = ac * mean_roes_j2(:,1);
mdls = ac *mean_roes_j2(:,2);
mdexs = ac *mean_roes_j2(:,3);
mdeys = ac *mean_roes_j2(:,4);
mdixs = ac *mean_roes_j2(:,5);
mdiys = ac *mean_roes_j2(:,6);

figure
subplot(1, 3, 1)
hold on
plot(dexs, deys)
plot(mdexs, mdeys)
xlabel('$a \delta e_x$ (m)', Interpreter='latex')
ylabel('$a \delta e_y$ (m)', Interpreter='latex')
axis square


subplot(1, 3, 2)
hold on
plot(dixs, diys, 'DisplayName','Osculating')
plot(mdixs, mdiys, 'DisplayName','Mean')
xlabel('$a \delta i_x$ (m)', Interpreter='latex')
ylabel('$a \delta i_y$ (m)', Interpreter='latex')
axis square
legend

subplot(1, 3, 3)
hold on
plot(dls, das)
plot(mdls, mdas)
xlabel('$a \delta \lambda$ (m)', Interpreter='latex')
ylabel('$a \delta a$ (m)', Interpreter='latex')

fontsize(11, 'points')

axis square