function [B] = chernick_control_matrix(singular_oe, mu)
    %{
    chernick_control_matrix calculates the Control Input Matrix defined in
    eqn 2.10, section 2.4
    
    Inputs:
        - singular_oe: keplerian orbit elements [a, e, i, RAAN, w, M] 
        - mu: Gravitational parameter of the primary attractor (m^3/s^2)
        
    Outputs: 
        - B: Control Input Matrix
    %}

% unpack
a = singular_oe(1);
e = singular_oe(2);
i = singular_oe(3);
% Omega = singular_oe(4);
w = singular_oe(5);
M = singular_oe(6);

nu = mean2true(M, e, 1e-6);

theta = w + nu;

n = mean_motion(mu, a);
eta = sqrt(1-e^2);
ex = e*cos(w);
ey = e*sin(w);
B = (1/(n*a)) * [
    2*e*sin(nu)/eta,            (2/eta)*(1+e*cos(nu)),                              0;
    -2*eta^2 / (1+e*cos(nu)),   0,                                                  0;
    eta*sin(theta),             eta*((2+e*cos(nu))*cos(theta) + ex)/(1+e*cos(nu)),  eta*ey*sin(theta)/(tan(i)*(1+e*cos(nu)));
    -eta*cos(theta),            eta*((2+e*cos(nu))*sin(theta) + ey)/(1+e*cos(nu)),  -eta*ex*sin(theta)/(tan(i)*(1+e*cos(nu)));
    0,                          0,                                                  eta*cos(theta)/(1+e*cos(nu));
    0,                          0,                                                  eta*sin(theta)/(1+e*cos(nu));
    ];

% u = nu + w;
% near-circular:
% B = (1/(n*a)) * [
%     0,          2,          0;
%     -2,         0,          0;
%     sin(u),     2*cos(u),   0;
%     -cos(u),    2*sin(u),   0;
%     0,          0,          cos(u);
%     0,          0,          sin(u);
% ];
end

