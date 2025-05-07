function Gamma = control_input_matrix(aoes, nu_k, mu)
%CONTROL_INPUT_MATRIX Computes Gamma
%   BEWARE: For use with MODIFIED QNS ROE!

a = aoes(1);
e = aoes(2);
i = aoes(3);
Omega = aoes(4);
w = aoes(5);
nu = aoes(6);

theta = nu_k + w;

n = mean_motion(mu, a);
eta = sqrt(1-e^2);
Gamma = (1/(n*a)) * [...
    2*e*sin(nu_k)/eta, 2*(1+e*cos(nu_k))/eta, 0;
    -2*eta^2/(1+e*cos(nu_k)),   0,    0;
    eta*sin(nu_k), eta* (e+cos(nu_k))*(2+e*cos(nu_k))/(1+e*cos(nu_k)), 0;
    -(eta/e) * cos(nu_k), (eta/e) * sin(nu_k) * (2+e*cos(nu_k)) / (1+e*cos(nu_k)), 0;
    0,  0,  eta*cos(theta)/(1+e*cos(nu_k));
    0,  0,  eta*sin(theta)/(1+e*cos(nu_k));
];

end