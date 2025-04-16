% OSC2MEAN converts a set of osculating Keplerian orbital elements to mean
% Keplerian orbital elements perturbed by Earth's oblateness using an
% iterative approach. Reference used is Spacecraft Formation Flying
% (Alfriend, 2010).
%
%   Inputs:
%     osc_elem - vector of osculating Keplerian orbital elements
%                a - semi-major axis [m]
%                e - eccentricity [-]
%                i - inclination [rad]
%                O - right ascension of the ascending node [rad]
%                w - argument of periapsis [rad]
%                M - mean anomaly [rad]
%     J2_flag - flag indicating whether or not to consider J2
%                1 -> J2 is enabled, calculate according to algorithm
%                0 -> J2 is disabled, osculating elements = mean elements
%                DEFAULT J2_flag = 1
% 
%   Outputs:
%     mean_elem - vector of mean Keplerian orbital elements
%                 a - semi-major axis [m]
%                 e - eccentricity [-]
%                 i - inclination [rad]
%                 O - right ascension of the ascending node [rad]
%                 w - argument of periapsis [rad]
%                 M - mean anomaly [rad]

function mean_elem = osc2mean(osc_elem,J2_flag)

    % Check inputs
    if (nargin < 2) || isempty(J2_flag)
        J2_flag = 1;
    end
    if (nargin < 1) || isempty(osc_elem)
        error('Must input mean elements set');
    end
    
    % Format input to column vector and set tolerance
    osc_elem = osc_elem(:);
    tol = 1e-08;
    
    % With J2, run iterative method
    if J2_flag == 1

        % Convert to osculating equinoctial elements
        osc_equioe = koe2equioe(osc_elem);
        
        % Convert to mean equinoctial elements
        equi_c_mean = osc2mean_NRiterator(osc_equioe, tol);
        
        % Convert to mean keplerian elements
        mean_elem = equioe2koe(equi_c_mean);
    
    % Without J2, elements are equal
    else
        
        mean_elem = osc_elem;
        
    end
    
end
