function dalpha_qns = quasi_nonsingular_roe(alpha_chief, alpha_deputy)
%QUASI_NONSINGULAR_ROE computes the quasi-nonsingular ROEs for a chief and deputy satellite.
% a: semi-major axis
% M: mean anomaly
% w: argument of periapsis
% e: eccentricity
% raan: right ascension of the ascending node
% i: inclination
% The c or d postfix in the variable name determines whether it is the
% chief or deputy spacecraft. 
ac = alpha_chief(1);
ec = alpha_chief(2);
ic = alpha_chief(3);
raanc = alpha_chief(4);
wc = alpha_chief(5);
nuc = alpha_chief(6);
Mc = true2mean(nuc, ec);

ad = alpha_deputy(1);
ed = alpha_deputy(2);
id = alpha_deputy(3);
raand = alpha_deputy(4);
wd = alpha_deputy(5);
nud = alpha_deputy(6);
Md = true2mean(nud, ed);

da = (ad-ac)/ac;
eta = sqrt(1-ec^2);
dlambda = Md - Mc + eta*(wd - wc + (raand - raanc)*cos(ic));
dex = ed*cos(wd) - ec*cos(wc);
dey = ed*sin(wd) - ec*sin(wc);
dix = id - ic;
diy = (raand - raanc)*sin(ic);

dalpha_qns = [da; dlambda; dex; dey; dix; diy;];
end