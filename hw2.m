%
close all
MU = 3.9860043550702260E+14; % m^3/s^2

function ds = ferm(~,s)
% Fundamental Equations of Relative Motion
% state vector s is comprized of:
x = s(1);
y = s(2);
z = s(3);
dx = s(4);
dy = s(5);
dz = s(6);
dtheta = s(7);
r0 = s(8);
dr0 = s(9);

MU = 3.9860043550702260E+14; % m^3/s^2
factor = -MU / (((r0+x)^2 + y^2 + z^2)^(3/2));

ddr0 = r0 * dtheta^2 - MU/r0^2;
ddtheta = -2*dr0*dtheta/r0;

ddx = factor * (r0 + x) + MU/(r0^2) + 2*dtheta*dy + ddtheta*y + dtheta^2*x;
ddy = factor * y - 2*dtheta*dx - ddtheta * x + dtheta^2 * y;
ddz = z * factor;

ds = [dx, dy,dz, ddx, ddy,ddz, ddtheta, dr0, ddr0].';
end


% initial Keplerian Orbital elements of Chief spacecraft (MMS-1):
a0 = 9.758384355377886E+07; % semi-major axis, meters
e0 = 8.820570082745063E-01; % eccentricity
i0 = 6.311884598781460E+01; % inclination, degrees
raan0 = 3.567045162674137E+02; % Right Ascencion of the Ascending Node (aka longitude of ascending node aka omega), degrees
argp0 = 1.582342516847609E+02; % Argument of periapsis (aka Argument of Perifocus), degrees
nu0 = 1.672699831765240E+02; % True Anomaly, degrees

da = 1E3;
da=0;
% initial Keplerian Orbital elements of Follower spacecraft
a1 = a0+da; % semi-major axis
e1 = e0-1e-3; % eccentricity
i1 =  i0 + 1E-3; % inclination
raan1 = raan0; %Right Ascencion of the Ascending Node
argp1 = argp0; % Argument of periapsis
nu1 = nu0-1e-2; % True Anomaly


dnu0 = (1+e0 * cos(deg2rad(nu0)))^2 * sqrt(MU / (a0^3 * (1-e0^2)^3));


% convert initial Keplerian orbital elements to cartesian
[r0, v0] = keplerian2ijk(a0, e0, i0, raan0, argp0, nu0);
[r1, v1] = keplerian2ijk(a1, e1, i1, raan1, argp1, nu1);

% compute relative initial conditions
rho_ECI = r1 - r0;
drho_ECI = v1 - v0;

% define intial RTN frame orientation in ECI
rhat = r0 / norm(r0);
h0 = cross(r0, v0);
hhat = h0 / norm(h0);
nhat = hhat;
that = cross(nhat, rhat);
rotation_matrix = [rhat, that, nhat];

% convert quantities into RTN
rho_RTN = rotation_matrix * rho_ECI;
omega_rtn_in_eci = nhat * dnu0;
drho_RTN = drho_ECI - cross(omega_rtn_in_eci,rho_ECI);
v0_RTN = rotation_matrix * v0;

% and compute intial conditions from there
x0 = rho_RTN(1);
y0 = rho_RTN(2);
z0 = rho_RTN(3);
vx0 = drho_RTN(1);
vy0 = drho_RTN(2);
vz0 = drho_RTN(3);
radius0 = norm(r0);
dradius0 = v0_RTN(1);
dtheta = dnu0;

initial_conditions = [x0, y0, z0, vx0, vy0, vz0, dtheta, radius0, dradius0].';


% propagate
num_orbits = 10;
Tp = 2*pi*sqrt(a0^3 / MU); % BEWARE: units must be seconds!

t_end = num_orbits*Tp; 
t = linspace(0, t_end, t_end / 30);

[~, s] = ode89(@ferm, t, initial_conditions);
rs = s(:,1);
ts = s(:,2);
ns = s(:,3);
vrs = s(:,4);
vts = s(:,5);
vns = s(:,6);




% now do it again but keplerian
nus0 = keplerian_propagation(nu0, e0, a0, t);
nus1 = keplerian_propagation(nu1, e1, a1, t);
states0 = zeros(6, length(t));
states1 = zeros(6, length(t));
for i = 1:length(t)
    [rijk0, vijk0] = keplerian2ijk(a0, e0, i0, raan0, argp0, rad2deg(nus0(i)));
    [rijk1, vijk1] = keplerian2ijk(a1, e1, i1, raan1, argp1, rad2deg(nus1(i)));
    state0 = [rijk0; vijk0];
    state1 = [rijk1; vijk1];

    states0(:, i) = state0;
    states1(:, i) = state1;
end

rs0 = states0(1:3, :);
rs1 = states1(1:3, :);
vs0 = states0(4:6, :);
vs1 = states1(4:6, :);
% calculate R, T, N vectors
hs0 = cross(rs0, vs0);
nhats0 = hs0./vecnorm(hs0);
rhats0 = rs0./vecnorm(rs0);
thats0 = cross(nhats0, rhats0);

% convert from ECI to RTN
rhos_ECI = rs1 - rs0;
drhos_ECI = vs1 - vs0;

rhos_RTN = zeros(size(rhos_ECI));
drhos_RTN = zeros(size(drhos_ECI));
for i = 1:length(t)
    rot_mat0 = [rhats0(:, i), thats0(:, i), nhats0(:, i)];
    rhos_RTN(:, i) = rot_mat0 * rhos_ECI(:, i);
    dnu = (1+e0 * cos(nus0(i)))^2 * sqrt(MU / (a0^3 * (1-e0^2)^3));
    omega_rtn_in_eci = nhats0(:,i) * dnu;
    drhos_RTN(:, i) =  drhos_ECI(:,i) - cross(omega_rtn_in_eci, rhos_ECI(:, i));
    
end

rs_kep = rhos_RTN(1,:);
ts_kep = rhos_RTN(2,:);
ns_kep = rhos_RTN(3,:);

vrs_kep = drhos_RTN(1,:);
vts_kep = drhos_RTN(2,:);
vns_kep = drhos_RTN(3,:);



figure;
plot_t = t / Tp;
subplot(4,2,1);
hold on
plot(rs / a0, ts/ a0);
plot(rs_kep/ a0, ts_kep / a0);
xlabel("x/a")
ylabel("y/a")
grid on;

subplot(4,2,3);
hold on
plot(ts / a0, ns/ a0);
plot(ts_kep / a0, ns_kep / a0);
xlabel("y/a")
ylabel("z/a")
grid on;


subplot(4,2,5);
hold on
plot(rs/ a0, ns/ a0);
plot(rs_kep / a0, ns_kep / a0);
xlabel("x/a")
ylabel("z/a")
grid on;


subplot(4,2,7);
hold on
plot(plot_t, rs, 'DisplayName', "R");
plot(plot_t, ts, 'DisplayName', "T");
plot(plot_t, ns, 'DisplayName', "N");
ylabel("Position [m]")
xlabel("Number of Orbits")
legend;
grid on;


plot_t = t / Tp;
subplot(4,2,2);
hold on
plot(vrs, vts);
plot(vrs_kep, vts_kep);
xlabel("Vx (m/s)")
ylabel("Vy (m/s)")
grid on;

subplot(4,2,4);
hold on
plot(vts, vns);
plot(vts_kep, vns_kep);
xlabel("Vy (m/s)")
ylabel("Vz (m/s)")
grid on;


subplot(4,2,6);
hold on
plot(vrs, vns);
plot(vrs_kep, vns_kep);
xlabel("Vx (m/s)")
ylabel("Vz (m/s)")
grid on;


subplot(4,2,8);
hold on
plot(plot_t, vrs, 'DisplayName', "Vr");
plot(plot_t, vts, 'DisplayName', "Vt");
plot(plot_t, vns, 'DisplayName', "Vn");
ylabel("Position [m]")
xlabel("Number of Orbits")
legend;
grid on;