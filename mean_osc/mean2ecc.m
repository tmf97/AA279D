% MEAN2ECC solves Kepler's equation for eccentric anomaly. Note that this
% function uses a Newton-Raphson method to numerically compute a value for
% eccentric anomaly.
%
%   Inputs:
%     M   - mean anomaly [rad]
%     e   - eccentricity [-]
%     tol - tolerance for Newton-Raphson iterator
%
%   Outputs:
%     E - eccentric anomaly [rad]

function E = mean2ecc(M, e, tol)

M = mod(M, 2*pi);

if M == 0 || M == pi
    % Return the known solutions (trivial)
    E = M;
else
    % Set up the problem based on an initial guess
    E0 = M;
    d = -(E0 - e*sin(E0) - M)/(1 - e*cos(E0));
    
    % Loop until the solution converges
    while abs(d) > tol
        E1 = E0 + d;
        
        d = -(E1 - e*sin(E1) - M)/(1 - e*cos(E1));
        
        E0 = E1;
    end
    
    E = E0;
end

end
