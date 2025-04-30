function pCi = pfCeci(asc, inc, argp)
%PFCECI Compute rotation matrix from perifocal frame to ECI

Casc = [cos(asc), -sin(asc), 0;
        sin(asc),  cos(asc), 0;
               0,         0, 1];

Cinc = [1,        0,        0;
        0, cos(inc), -sin(inc);
        0, sin(inc), cos(inc)];

Cargp = [cos(argp), -sin(argp), 0; 
         sin(argp),  cos(argp), 0; 
                 0,          0, 1];

pCi = Casc * Cinc * Cargp;
end