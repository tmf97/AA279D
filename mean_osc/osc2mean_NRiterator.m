% OSC2MEAN_NRITERATOR iteratively computes a set of mean equinoctial orbital
% elements from osculating orbital elements using an NR iteration containing the
% Jacobian of the osculating-to-mean transformation.
% 
%   Inputs:
%     osc_equi_elem - vector of osculating equinoctial orbital elements
%                     a   - semi-major axis [m]
%                     Psi - mean longitude [rad]
%                     tq1 - [-]
%                     tq2 - [-]
%                     p1  - [-]
%                     p2  - [-]
%     tol           - iterative solver tolerance (typically 1e-12)
% 
%   Outputs:
%     mean_equi_elem - vector of mean equinoctial orbital elements
%                      a   - semi-major axis [m]
%                      Psi - mean longitude [rad]
%                      tq1 - [-]
%                      tq2 - [-]
%                      p1  - [-]
%                      p2  - [-]

function mean_equi_elem = osc2mean_NRiterator(osc_equi_elem, tol)

mean_equi_elem = osc_equi_elem;
R = 1;
niter = 0;
while abs(R) > tol
  niter = niter+1;
  [~, osc_loop, ~] = transformationmatrix_osc2mean_equinoctial(mean_equi_elem);
  delta = osc_equi_elem - osc_loop;
  R = norm(delta, inf);
  mean_equi_elem = mean_equi_elem + delta;
  if niter > 100
      disp('Osc2Mean iterations > 100');
      break;
  end
  
end

end
