function AOE_deputy = roe2aoe(AOE_chief, ROE_deputy)
%ROE2AOE Relative to Absolute Orbital Elements
%   Computes the absolute orbital elements of the deputy spacecraft from
%   Keplierian absolute oribital elements of the chief, and
%   quasi-nonsingular relative orbital elements of the deputy.

% unpack chief array
ac = AOE_chief(1);
ec = AOE_chief(2);
ic = AOE_chief(3);
Omegac = AOE_chief(4); % RAAN
wc = AOE_chief(5); % argp
nuc = AOE_chief(6); % true anomaly
uc = wc + nuc; % argument of latitude of chief

% unpack deputy array
% angular quantities are re-normalized by ac
da = ROE_deputy(1);
dlambda = ROE_deputy(2) / ac;
dex = ROE_deputy(3) / ac;
dey = ROE_deputy(4) / ac;
dix = ROE_deputy(5) / ac;
diy = ROE_deputy(6) / ac;

% compute deputy absolute orbital elements
ad = ac + da;
id = ic + dix;
Omegad = Omegac + diy / sin(ic);

wd = atan2(dey + ec*sin(wc), dex + ec*cos(wc));

% want to avoid numerical error when sin(wd) is very small
if abs(sin(wd)) < 0.1
    ed = (dex + ec*cos(wc)) / cos(wd);
else
    ed = (dey + ec*sin(wc)) / sin(wd);
end

delta_u = dlambda - (Omegad - Omegac) * cos(ic);
ud = uc + delta_u;
nud = ud - wd;

AOE_deputy = [ad, ed, id, Omegad, wd, nud];
end