close all
clear truth_ephem
clear prop_ephem
% validate propagator:
MU = 3.9860043550702260E+14; % m^3/s^2
opts = odeset("RelTol",1e-13,"AbsTol",1e-15);

% read truth ephemerides:
% one column vector for all state.
[mms1_ephem, ts] = read_horizons('mms-1.txt');
[mms2_ephem, ~]  = read_horizons('mms-2.txt');
[mms3_ephem, ~]  = read_horizons('mms-3.txt');
[mms4_ephem, ~]  = read_horizons('mms-4.txt');

start_n = 500;
mms1_ephem = mms1_ephem(:, start_n:end);
mms2_ephem = mms2_ephem(:, start_n:end);
mms3_ephem = mms3_ephem(:, start_n:end);
mms4_ephem = mms4_ephem(:, start_n:end);

ts = ts(start_n:end);

truth_ephem(:,:,1) = mms1_ephem; 
truth_ephem(:,:,2) = mms2_ephem;
truth_ephem(:,:,3) = mms3_ephem;
truth_ephem(:,:,4) = mms4_ephem;



% convert julian days to delta seconds
enhancement = 2;
ts = (ts - ts(1)) * 24 * 3600;
prop_ts = ts(end) * linspace(0, 1, enhancement*length(ts) - (enhancement - 1));
figure
for i=1:4
    ephem = truth_ephem(:,:,i);
    koe_ic = ephem(:, 1);
    [ri, vi] = koe2pv(koe_ic, MU);
    init_state = [ri; vi];
    [~, propagated_states] = ode45(@newton3d_J2_moon, prop_ts, init_state, opts);
    propagated_states = propagated_states.';
    truth_states = zeros(size(ephem));
    for j=1:length(ephem)
        [rt, vt] = koe2pv(ephem(:, j), MU);
        truth_states(:, j) = [rt; vt];
    end
    
    if false
        subplot(4, 3, 3*i-2);
        hold on
        plot(ts, truth_states(1,:), 'DisplayName', "MMS-" + i + " Truth");
        plot(prop_ts, propagated_states(1,:), 'DisplayName', "MMS-" + i + " Propagated");
        hold off
        ylabel("[m]")

        subplot(4, 3, 3*i-1);
        hold on
        plot(ts, truth_states(2,:), 'DisplayName', "MMS-" + i + " Truth");
        plot(prop_ts, propagated_states(2,:), 'DisplayName', "MMS-" + i + " Propagated");
        hold off
        ylabel("[m]")


        subplot(4, 3, 3*i);
        hold on
        plot(ts, truth_states(3,:), 'DisplayName', "MMS-" + i + " Truth");
        plot(prop_ts, propagated_states(3,:), 'DisplayName', "MMS-" + i + " Propagated");
        hold off
        ylabel("[m]")        

    else
        subplot(4, 2, -1+2*i);
        prop_pos = propagated_states(1:3, 1:enhancement:end);
        prop_vel = propagated_states(4:6, 1:enhancement:end);
        pos_errors = prop_pos - truth_states(1:3, :);
        vel_errors = prop_vel - truth_states(4:6, :);
        plot(ts, vecnorm(pos_errors));

        ylabel("MMS-" + i + " Error [m]")
        yscale('log')
        subplot(4, 2, 2*i);
        plot(ts, vecnorm(vel_errors));
        ylabel("MMS-" + i + " Error [m/s]")
        yscale('log')

    end

    prop_ephem(:,:,i) = propagated_states;
end
xlabel("Time")


mms2_rtns = zeros(3, length(ts));
mms3_rtns = zeros(3, length(ts));
mms4_rtns = zeros(3, length(ts));

a0 = truth_ephem(1,1,1);

for i = 1:length(prop_ts)
    [mms2_rtns(:,i), ~] = pvs2rtn(prop_ephem(1:3,i,1), prop_ephem(4:6,i,1), prop_ephem(1:3,i,2), prop_ephem(4:6,i,2), MU, a0);
    [mms3_rtns(:,i), ~] = pvs2rtn(prop_ephem(1:3,i,1), prop_ephem(4:6,i,1), prop_ephem(1:3,i,3), prop_ephem(4:6,i,3), MU, a0);
    [mms4_rtns(:,i), ~] = pvs2rtn(prop_ephem(1:3,i,1), prop_ephem(4:6,i,1), prop_ephem(1:3,i,4), prop_ephem(4:6,i,4), MU, a0);
end

Rs2 = mms2_rtns(1, :);
Ts2 = mms2_rtns(2, :);
Ns2 = mms2_rtns(3, :);

Rs3 = mms3_rtns(1, :);
Ts3 = mms3_rtns(2, :);
Ns3 = mms3_rtns(3, :);

Rs4 = mms4_rtns(1, :);
Ts4 = mms4_rtns(2, :);
Ns4 = mms4_rtns(3, :);



figure
subplot(3,2,1);
hold on
plot(prop_ts, Rs2);
plot(prop_ts, Rs3);
plot(prop_ts, Rs4);
ylabel("R (m)")
grid on;


subplot(3,2,3);
hold on
plot(prop_ts, Ts2);
plot(prop_ts, Ts3);
plot(prop_ts, Ts4);
ylabel("T (m)")
grid on;


subplot(3,2,5);
hold on
plot(prop_ts, Ns2);
plot(prop_ts, Ns3);
plot(prop_ts, Ns4);
ylabel("N (m)")
xlabel("Orbital Periods")
grid on;


subplot(2,4,3)
hold on
plot(Ts2, Rs2)
plot(Ts3, Rs3)
plot(Ts4, Rs4)
xlabel("T (m)")
ylabel("R (m)")
grid on;
axis square;

subplot(2,4,4)
hold on
plot(Ns2, Rs2)
plot(Ns3, Rs3)
plot(Ns4, Rs4)
xlabel("N (m)")
ylabel("R (m)")
grid on;
axis square;

subplot(2,4,7)
hold on
plot(Ts2, Ns2)
plot(Ts3, Ns3)
plot(Ts4, Ns4)
xlabel("T (m)")
ylabel("N (m)")
grid on;
axis square;

subplot(2,4,8)
hold on
plot3(Rs2, Ts2, Ns2)
plot3(Rs3, Ts3, Ns3)
plot3(Rs4, Ts4, Ns4)
xlabel("R (m)")
ylabel("T (m)")
zlabel("N (m)")
grid on;
axis square;

%%%%%
% The TRUTH
%%%%%

tmms2_rtns = zeros(3, length(ts));
tmms3_rtns = zeros(3, length(ts));
tmms4_rtns = zeros(3, length(ts));

a0 = truth_ephem(1,1,1);

for i = 1:length(ts)
    tmms1_koe = mms1_ephem(:,i);
    tmms2_koe = mms2_ephem(:,i);
    tmms3_koe = mms3_ephem(:,i);
    tmms4_koe = mms4_ephem(:,i);
    [tmms1_r, tmms1_v] = koe2pv(tmms1_koe, MU);
    [tmms2_r, tmms2_v] = koe2pv(tmms2_koe, MU);
    [tmms3_r, tmms3_v] = koe2pv(tmms3_koe, MU);
    [tmms4_r, tmms4_v] = koe2pv(tmms4_koe, MU);

    [tmms2_rtns(:,i), ~] = pvs2rtn(tmms1_r, tmms1_v, tmms2_r, tmms2_v, MU, a0);
    [tmms3_rtns(:,i), ~] = pvs2rtn(tmms1_r, tmms1_v, tmms3_r, tmms3_v, MU, a0);
    [tmms4_rtns(:,i), ~] = pvs2rtn(tmms1_r, tmms1_v, tmms4_r, tmms4_v, MU, a0);
end

Rs2 = tmms2_rtns(1, :);
Ts2 = tmms2_rtns(2, :);
Ns2 = tmms2_rtns(3, :);

Rs3 = tmms3_rtns(1, :);
Ts3 = tmms3_rtns(2, :);
Ns3 = tmms3_rtns(3, :);

Rs4 = tmms4_rtns(1, :);
Ts4 = tmms4_rtns(2, :);
Ns4 = tmms4_rtns(3, :);



figure
title('Truth')
subplot(3,2,1);
hold on
plot(ts, Rs2);
plot(ts, Rs3);
plot(ts, Rs4);
ylabel("R (m)")
grid on;


subplot(3,2,3);
hold on
plot(ts, Ts2);
plot(ts, Ts3);
plot(ts, Ts4);
ylabel("T (m)")
grid on;


subplot(3,2,5);
hold on
plot(ts, Ns2);
plot(ts, Ns3);
plot(ts, Ns4);
ylabel("N (m)")
xlabel("Orbital Periods")
grid on;


subplot(2,4,3)
hold on
plot(Ts2, Rs2)
plot(Ts3, Rs3)
plot(Ts4, Rs4)
xlabel("T (m)")
ylabel("R (m)")
grid on;
axis square;

subplot(2,4,4)
hold on
plot(Ns2, Rs2)
plot(Ns3, Rs3)
plot(Ns4, Rs4)
xlabel("N (m)")
ylabel("R (m)")
grid on;
axis square;

subplot(2,4,7)
hold on
plot(Ts2, Ns2)
plot(Ts3, Ns3)
plot(Ts4, Ns4)
xlabel("T (m)")
ylabel("N (m)")
grid on;
axis square;

subplot(2,4,8)
hold on
plot3(Rs2, Ts2, Ns2)
plot3(Rs3, Ts3, Ns3)
plot3(Rs4, Ts4, Ns4)
xlabel("R (m)")
ylabel("T (m)")
zlabel("N (m)")
grid on;
axis square;


figure
grid on
hold on
axis equal
Re = 6378000; % m
[X,Y,Z] = sphere;
surf(X*Re, Y*Re, Z*Re, 'DisplayName', 'Earth')
plot3(truth_states(1,:), truth_states(2,:), truth_states(3,:))
plot3(propagated_states(1,1:enhancement:end), propagated_states(2,1:enhancement:end), propagated_states(3,1:enhancement:end))
plot3(propagated_states(1,:), propagated_states(2,:), propagated_states(3,:))


figure 
hold on
plot(ts, vecnorm(truth_states(1:3,:)));
plot(prop_ts, vecnorm(propagated_states(1:3,:)));
ylabel("[m]")

hs_true = zeros(size(ts));
for i=1:length(ts)
    r = truth_states(1:3,i);
    v = truth_states(4:6,i);
    h = cross(r,v);
    hs_true(i) = norm(h);
end

hs_prop = zeros(size(prop_ts));
for i=1:length(prop_ts)
    rp = propagated_states(1:3,i);
    vp = propagated_states(4:6,i);
    hp = cross(rp, vp);
    hs_prop(i) = norm(hp);
end

figure 
hold on
plot(ts, hs_true);
plot(prop_ts, hs_prop);
ylabel("h")


