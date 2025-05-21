function B = control_input_matrix_v2(aoes, mu)
%CONTROL_INPUT_MATRIX_V2 Computes the control input matrix
%   From equation 4 of Steindorf 2017. Uses AOEs of the chief
%   It's v2 because i forgot i implemented it already
a = aoes(1);
e = aoes(2);
i = aoes(3);
% Omega = aoes(4);
w = aoes(5);
f = aoes(6);

n = mean_motion(mu, a);
eta = sqrt(1-e^2);
ex = e*cos(w);
ey = e*sin(w);

B = (1/(a*n)) * [
  2*(1 + e*cos(f))/eta,                             0;
  eta*((2 + e*cos(f))*cos(w+f) + ex)/(1+e*cos(f)),  eta*ey*sin(w+f)/(tan(i)*(1+e*cos(f)))
  eta*((2 + e*cos(f))*sin(w+f) + ey)/(1+e*cos(f)),  -eta*ex*sin(w+f)/(tan(i)*(1+e*cos(f)))
  0,                                                eta*cos(w+f) / (1+e*cos(f));
  0,                                                eta*sin(w+f) / (1+e*cos(f));
];
end