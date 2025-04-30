function ds = newton3d(~,s)
%NEWTON3D Orbital propagator based off Newtonian gravity propagation
%   s is the state vector of position and velocity in meters and m/s
MU = 3.9860043550702260E+14; % m^3/s^2

r = s(1:3);
v = s(4:6);
ddr = -MU * r/(norm(r)^3);
ds = [v; ddr];
end
