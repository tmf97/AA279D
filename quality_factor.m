function Q = quality_factor(sep12, sep13, sep14, target_sep)
%QUALITY_FACTOR Computes MMS quality factor
%   Calculations described in Fuselier, S.A., Lewis, W.S., Schiff, C. et
%   al. Magnetospheric Multiscale Science Mission Profile and Operations.
%   Space Sci Rev 199, 77–103 (2016).
%   https://doi.org/10.1007/s11214-014-0087-x

mag_sep12 = norm(sep12);
mag_sep13 = norm(sep13);
mag_sep14 = norm(sep14);
L = mean([mag_sep12, mag_sep13, mag_sep14]);
real_volume = tetrahedron_volume(sep12, sep13, sep14);

% volume of a regular tetrahedron
ideal_volume = L^3 / (6*sqrt(2));

Qv = real_volume / ideal_volume;

lookup_table = dictionary( ...
    10, [4, 6, 18, 24], ...
    25, [15, 20, 35, 40], ...
    30, [19.313, 23.15, 42.075, 49.475], ...
    40, [25, 30, 55, 65], ...
    60, [45, 50, 75, 80], ...
    160, [135, 140, 190, 210], ...
    400, [250, 300, 550, 600]);

Ls = lookup_table(target_sep);
L1 = Ls(1);
L2 = Ls(2);
L3 = Ls(3);
L4 = Ls(4);

if L < L1 || L > L4
    Qs = 0;
end
if L2 < L && L < L3
    Qs = 1;
end
if L1 < L && L < L2
    Qs = (L-L1)^2 * (L+L1-2*L2)^2/(L2 - L1)^4;
end
if L3 < L && L < L4
    Qs = (L-L4)^2 * (L-2L3+L4)^2 / (L2-L1)^4;
end

Q = Qs * Qv;
end