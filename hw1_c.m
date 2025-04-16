close all

% section C

% 3D 2 Body Problem
MU = 3.9860043550702260E+05; % km^3/s^2
r0 = [125365.4; -47010.55; -78368.37]; % km
v0 = [1.0252; -0.04922; 0.01911]; % km/s
s0 = [r0; v0];

num_orbits = 2;
t_end = num_orbits*24*60*60;
t = linspace(0, t_end, t_end / 2);

[~, s] = ode45(@newton3d, t, s0);
[~, s_j2] = ode45(@newton3d_J2, t, s0);
figure;
plot_t = t / (24* 3600);
subplot(3,2,1);
hold on
plot(plot_t, s(:, 1), "DisplayName", "No Perturbation")
plot(plot_t, s_j2(:, 1), "DisplayName", "J2 Perturbations")
ylabel("X Position [km]")
legend;
grid on;

subplot(3,2,3);
hold on
plot(plot_t, s(:, 2), "DisplayName", "No Perturbation")
plot(plot_t, s_j2(:, 2), "DisplayName", "J2 Perturbations")
ylabel("Y Position [km]")
legend;
grid on;

subplot(3,2,5);
hold on
plot(plot_t, s(:, 3), "DisplayName", "No Perturbation")
plot(plot_t, s_j2(:, 3), "DisplayName", "J2 Perturbations")
ylabel("Z Position [km]")
xlabel("Time [days]")
legend;
grid on;


subplot(3,2,2);
hold on
plot(plot_t, s(:, 4), "DisplayName", "No Perturbation")
plot(plot_t, s_j2(:, 4), "DisplayName", "J2 Perturbations")
ylabel("X Velocity [km/s]")
legend;
grid on;

subplot(3,2,4);
hold on
plot(plot_t, s(:, 5), "DisplayName", "No Perturbation")
plot(plot_t, s_j2(:, 5), "DisplayName", "J2 Perturbations")
ylabel("Y Velocity [km/s]")
legend;
grid on;

subplot(3,2,6);
hold on
plot(plot_t, s(:, 6), "DisplayName", "No Perturbation")
plot(plot_t, s_j2(:, 6), "DisplayName", "J2 Perturbations")
ylabel("Z Velocity [km/s]")
xlabel("Time [days]")
legend;
grid on;

% 3D plot of orbital positions
figure
hold on
Re = 6378.1; % km
[X,Y,Z] = sphere;
surf(X*Re, Y*Re, Z*Re, 'DisplayName', 'Earth')
plot3(s(:, 1), s(:,2), s(:,3), 'DisplayName', "Without J2")
plot3(s_j2(:, 1), s_j2(:,2), s_j2(:,3), 'DisplayName', "With J2")
legend
grid on

axis equal
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('MMS-1 Orbit Propagation via Newton Gravity Propagation')
hold off




% s is the state vector of position and velocity
function ds = newton3d(~,s)
MU = 3.9860043E+05; % km^3/s^2
r = s(1:3);
v = s(4:6);
ddr = -MU * r/(norm(r)^3);
ds = [v; ddr];
end


function ds = newton3d_J2(~,s)
Re = 6378.1; % km
MU = 3.9860043550702260E+05; % km^3/s^2
r = s(1:3);
x = r(1);
y = r(2);
z = r(3);
r_norm = norm(r);
v = s(4:6);
a_gravity = -MU * r/(r_norm^3);

J2 = 0.00108263;


factor = (3/2) * J2*MU*Re^2 / (r_norm^5);

common = (5 * z^2 / r_norm^2) - 1;
a_J2 = factor *[x * common; y * common; z * (common - 2)];

a_tot = a_J2 + a_gravity;
ds = [v; a_tot];
end



% section D
function En = mean_to_eccentric_anomaly(M, e, tol)
% Use Newton-Raphson iteration to solve Kepler's equation E - e*sin(E) = M
% for the eccentric anomaly (E) given the eccentricity (e) and the mean
% anomaly (M).
En = M + e * sin(M); % initial guess for E
while abs(M - (En - e*sin(En))) > tol
    Enp1 = En - ((M-En + e*sin(En))./(e*cos(En) - 1));
    En=Enp1;
end
end

function M = eccentric_to_mean_anomaly(E, e)
M = E - e * sin(E);
end

function nu = eccentric_to_true_anomaly(E, e)
% use Kepler's equations to solve for True Anomaly (nu) from the Eccentric
% Anomaly (E) and the eccentricity (e).

nu = mod(atan2(sqrt(1-e^2)*sin(E), cos(E)-e), 2*pi);

end

function E = true_to_eccentric_anomaly(nu, e)
% use Kepler's equations to solve for Eccentric Anomaly (E) from the True
% Anomaly (nu) and the eccentricity (e).

E = atan2(sin(nu)*sqrt(1-e^2), cos(nu) + e);

end

% how _do_ you write unittests in matlab??
E_test = 0.01; % rad
e_test = 0.1;
nu_test = eccentric_to_true_anomaly(E_test, e_test);
E_test_2 = true_to_eccentric_anomaly(nu_test, e_test);
assert((E_test - E_test_2)<1e-7)




MU = 3.9860043550702260E+14; % m^3/s^2


DEG_2_RAD = pi / 180;
RAD_2_DEG = 180 / pi;
% take initial keplerian orbit
% convert cartesian into keplerian

[a, e, inc, RAAN, omega, nu] = ijk2keplerian(r0*1000, v0*1000);
M = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu * DEG_2_RAD, e),e);

% propagate orbit
n = mean_motion(MU, a);
Ms = n.*t + M;
Es = mean_to_eccentric_anomaly(Ms, e, 1e-12);
nus = eccentric_to_true_anomaly(Es, e);

states = zeros(6, length(t));
for i = 1:length(t)
[r_ijk, v_ijk] = keplerian2ijk(a, e, inc, RAAN, omega, nus(i) * RAD_2_DEG);
s_ijk = [r_ijk; v_ijk];
states(:, i) = s_ijk;
end

states = states / 1000;

errors = s.' - states;
errors = errors.';

figure;
plot_t = t / (24* 3600);
subplot(3,2,1);
hold on
plot(plot_t, abs(errors(:, 1)))
ylabel("X Error [m]")
yscale log
grid on;

subplot(3,2,3);
hold on
plot(plot_t, abs(errors(:, 2)))
ylabel("Y Error [m]")
yscale log
grid on;

subplot(3,2,5);
hold on
plot(plot_t, abs(errors(:, 3)))
ylabel("Z Error [m]")
xlabel("Time [days]")
yscale log
grid on;


subplot(3,2,2);
hold on
plot(plot_t, abs(errors(:, 4)))
ylabel("XV Error [m/s]")
yscale log
grid on;

subplot(3,2,4);
hold on
plot(plot_t, abs(errors(:, 5)))
ylabel("YV Error [m/s]")
yscale log
grid on;

subplot(3,2,6);
hold on
plot(plot_t, abs(errors(:, 6)))
ylabel("ZV Error [m/s]")
xlabel("Time [days]")
yscale log
legend;
grid on;


% 3D plot of orbital positions
figure
hold on
Re = 6378.1; % km
[X,Y,Z] = sphere;
surf(X*Re, Y*Re, Z*Re, 'DisplayName', 'Earth')
plot3(s(:, 1), s(:,2), s(:,3), 'DisplayName', "Newton")
plot3(states(1,:), states(2,:), states(3,:), 'DisplayName', "Keplerian")
legend
grid on

axis equal
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('MMS-1 Orbit Propagation via Newton Gravity Propagation')
hold off