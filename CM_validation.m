% starting a new file so I don't screw myself over with messing up a file
% that does work the day before the final is due. Bad coding practice, I
% know, but I'm desperate.

close all
% Constants
J2 = 0.0010826358191967; % Earth's J2 coefficient
mu = 3.986004415e14; % Earth gravitational parameter [m^3/s^2]
rE = 6.378136300e6; % Earth radius [m]

% chief AOEs
AOE_chief = [63781000, 0.6, deg2rad(50.11885), deg2rad(56.7045), deg2rad(58.2343), deg2rad(12)];
ROEs = [20, 2000, 50, 100, 30, 200] / AOE_chief(1); % meters

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

ts = linspace(0, 10*Tp, 1e4);
ts_2 = linspace(ts(end), ts(end)+10*Tp, 1e4);

dt = ts(2) - ts(1);

opts = odeset("RelTol",5e-14,"AbsTol",1e-17);


[~, states_chief_j2] = ode45(@newton3d_J2, ts, sc, opts);
[~, states_deputy_j2] = ode45(@newton3d_J2, ts, sd, opts);

s_chief = states_chief_j2(end, :);
s_deputy = states_deputy_j2(end, :);
% apply maneuver to deputy. Positive cross-track burn should increase SMA
maneuver_dv_rtn = [0; 0; 0.01];

maneuver_dv_eci = rtn2eci(s_chief(1:3), s_chief(4:6), maneuver_dv_rtn);
s_deputy_maneuvered = s_deputy.' + [0; 0; 0; maneuver_dv_eci];

% propagate post-maneuver
[~, states_chief_j2_post_maneuver] = ode45(@newton3d_J2, ts_2, s_chief, opts);
[~, states_deputy_j2_post_maneuver] = ode45(@newton3d_J2, ts_2, s_deputy_maneuvered, opts);

all_states_chief = [states_chief_j2; states_chief_j2_post_maneuver];
all_states_deputy = [states_deputy_j2; states_deputy_j2_post_maneuver];

kep_states_chief_j2 = zeros(size(all_states_chief));
kep_states_deputy_j2 = zeros(size(all_states_chief));
roes_j2 = zeros(size(all_states_chief));
mean_oes_chief_j2 = zeros(size(all_states_chief));
mean_oes_deputy_j2 = zeros(size(all_states_chief));
mean_roes_j2 = zeros(size(all_states_chief));
rtns_j2 = zeros(size(all_states_chief));

for i=1:length(all_states_chief)
    state_chief_j2= all_states_chief(i,:);
    state_deputy_j2 = all_states_deputy(i,:);
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
    [rtn, rtn_dot] = pvs2rtn(r0, v0,r1, v1);
    rtns_j2(i,:) = [rtn; rtn_dot];

end


% calculate STM propagation
stm_roes = zeros(size(all_states_deputy));
stm_roes(1,:) = ROEs.';

for i = 2:length(ts)
    % STM calcs
    stm = chernick_J2_stm(kep_states_chief_j2(i,:), dt, rE, MU, J2);
    stm_roe = stm * stm_roes(i-1,:).';
    stm_roes(i,:) = stm_roe;
end

% apply maneuver via control input matrix

phi = chernick_J2_stm(kep_states_chief_j2(i,:), dt, rE, MU, J2);
koe = kep_states_chief_j2(i,:);

e = koe(2);
true_anom = koe(6);
mean_anom = true2mean(true_anom, e);
B = chernick_control_matrix([koe(1:5), mean_anom], MU);

stm_roes(i,:)
roe_post_maneuver = phi * stm_roes(i,:).' + B * maneuver_dv_rtn

stm_roes(1+length(ts),:) = roe_post_maneuver;

for i = length(ts)+2:length(ts)+length(ts_2)
    % STM calcs
    stm = chernick_J2_stm(kep_states_chief_j2(i,:), dt, rE, MU, J2);
    stm_roe = stm * stm_roes(i-1,:).';
    stm_roes(i,:) = stm_roe;
end



% 
% as = kep_states_chief_j2(:,1);
% es = kep_states_chief_j2(:,2);
% is = kep_states_chief_j2(:,3);
% Omegas = kep_states_chief_j2(:,4);
% ws = kep_states_chief_j2(:,5);
% nus = kep_states_chief_j2(:,6);
% 
% mean_as = mean_oes_chief_j2(:,1);
% mean_es = mean_oes_chief_j2(:,2);
% mean_is = mean_oes_chief_j2(:,3);
% mean_Omegas = mean_oes_chief_j2(:,4);
% mean_ws = mean_oes_chief_j2(:,5);
% mean_nus = mean_oes_chief_j2(:,6);

roes_j2 = roes_j2 * ac;
mean_roes_j2 = mean_roes_j2 * ac;

das = roes_j2(:,1);
dlambdas = roes_j2(:,2);
dexs = roes_j2(:,3);
deys = roes_j2(:,4);
dixs = roes_j2(:,5);
diys = roes_j2(:,6);

mean_das = mean_roes_j2(:,1);
mean_dlambdas = mean_roes_j2(:,2);
mean_dexs = mean_roes_j2(:,3);
mean_deys = mean_roes_j2(:,4);
mean_dixs = mean_roes_j2(:,5);
mean_diys = mean_roes_j2(:,6);


% denormalize STM ROEs
stm_roes = stm_roes * ac;

figure
subplot(2, 3, 1)
hold on
plot(dexs, deys)
% plot(mean_dexs, mean_deys)
plot(stm_roes(:,3), stm_roes(:,4))
xlabel('$a \delta e_x$ (m)', Interpreter='latex')
ylabel('$a \delta e_y$ (m)', Interpreter='latex')
axis square


subplot(2, 3, 2)
hold on
plot(dixs, diys, 'DisplayName','Osculating')
% plot(mean_dixs, mean_diys, 'DisplayName','Mean')
plot(stm_roes(:,5), stm_roes(:,6), 'DisplayName','STM')
xlabel('$a \delta i_x$ (m)', Interpreter='latex')
ylabel('$a \delta i_y$ (m)', Interpreter='latex')
axis square
legend

subplot(2, 3, 3)
hold on
plot(dlambdas, das)
% plot(mean_dlambdas, mean_das)
plot(stm_roes(:,2), stm_roes(:,1))
xlabel('$a \delta \lambda$ (m)', Interpreter='latex')
ylabel('$a \delta a$ (m)', Interpreter='latex')
axis square



subplot(2, 3, 4)
hold on
plot(stm_roes(:,3) - dexs, stm_roes(:,4) - deys)
xlabel('$a \delta e_x$ Error (m)', Interpreter='latex')
ylabel('$a \delta e_y$ Error (m)', Interpreter='latex')
axis square equal
pbaspect([1,1,1])

subplot(2, 3, 5)
hold on
plot(stm_roes(:,5)-dixs, stm_roes(:,6)-diys)
xlabel('$a \delta i_x$ Error (m)', Interpreter='latex')
ylabel('$a \delta i_y$ Error (m)', Interpreter='latex')
axis square equal
pbaspect([1,1,1])


subplot(2, 3, 6)
hold on
plot(stm_roes(:,2) - dlambdas, stm_roes(:,1) - das)
xlabel('$a \delta \lambda$ Error (m)', Interpreter='latex')
ylabel('$a \delta a$ Error (m)', Interpreter='latex')
axis square equal
pbaspect([1,1,1])

fontsize(14, 'points')

