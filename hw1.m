function [r,v] = kepler2cartesian(a,e,I,w,O,mu,t,tp)    
    n = sqrt(mu/a^3); % mean motion
    % asssume t and tp are in julian days
    dt = (t-tp) * 24*60*60;
    M = n*dt; % mean anomaly
    M = mod(M, 2*pi);
    E = netwon_raphson(M,e,eps(2*pi)); % eccentric anomaly

    % ta = 2 * atan(sqrt((1+e)/(1-e)) * tan(E/2)); % true anomaly
    % r = a * ((1-e^2)/(1+e*cos(ta))); % distance
    
    % calculate pos/vel in perifocal frame using Kepler's equations
    xp = a*(cos(E)-e);
    b = a*sqrt(1-e^2);
    yp = b * sin(E);
    zp = 0;
    
    Edot = n ./ (1-e*cos(E));   

    vxp = -a * Edot.*sin(E);
    vyp = b*Edot.*cos(E);
    vzp = 0;

    rp = [xp; yp; zp];
    vp = [vxp; vyp; vzp];
    pci = pfCeci(O, I, w);
    r = pci * rp;
    v = pci * vp;

end


    % transform from perifocal frame to ECI
    function pCi = pfCeci(asc, inc, argp)
       Casc = [cos(asc), sin(asc), 0;
               -sin(asc), cos(asc), 0;
               0,0,1];
       Cinc = [1,0,0;
               0, cos(inc), sin(inc);
               0, -sin(inc), cos(inc)];
       Cargp = [cos(argp), sin(argp), 0; 
                -sin(argp), cos(argp), 0; 
                0,0,1];

       pCi = Casc * Cinc * Cargp;
    end

    function En = netwon_raphson(m, e, tol)
        En = m + e * sin(m); % initial guess for E
        while abs(m - (En - e*sin(En))) > tol
            Enp1 = En - ((m-En + e*sin(En))./(e*cos(En) - 1));
            En=Enp1;
        end
    end



    
e = 8.820570082745063E-01;
lon = 3.567045162674137E+02 * pi / 180; % rad
sma = 9.758384355377886E+04; % km
inc = 6.311884598781460E+01 * pi /180; 
argp = 1.582342516847609E+02 *pi / 180;
ta = 1.672699831765240E+02 * pi /180;

t = 2460766.500000000;
tp = 2460765.579404851887;
MU = 3.9860043550702260E+05; % km^3/s^2

[r,v] = kepler2cartesian(sma, e, inc, argp, lon, MU, t, tp);
disp(r);
disp(v);