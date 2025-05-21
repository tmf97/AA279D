function [t_maneuvers, manuevers, total_cost] = impulsive_control(chief_oe, initial_aroe, desired_aroe, dt, STM, CM, Rp, mu, J2)
    %{
    impulsive_control calculates the optimal maneuver times and magnitudes
    to achieve the given difference in ROE states for the provided 
    reconfiguration time. This code is based on Chernick's PhD Thesis. It
    assumes a near-circular chief orbit. 

    Generally, we follow Algorithm 6.1 applied in the near-circular case. 

    This function is designed around SI units (m, s, radians). Use other 
    units at your own risk. 
    
    Inputs:
        - chief_oe: chief's keplerian orbit elements [a, e, i, RAAN, w, M] 
        - initial_aroe: deputy's initial dimensionalized QNS ROE state 
                a_c[da, dlambda, dex, dey, dix, diy]
        - desired_aroe: desired dimensionalized QNS ROE state
                a_c[da, dlambda, dex, dey, dix, diy]
        - dt: reconfiguration window time (s)
        - STM: function which given the chief OE's, and a time will return
                the STM by which to multiply ROEs to get the current ROE 
                state
        - CM: returns the control input matrix given the chief's ROEs
        - Rp: radius of the primary attractor (m)
        - mu: Gravitational parameter of the primary attractor (m^3/s^2)
        - J2: J2 coefficient of the primary attractor
        
    Outputs: 
        - t_maneuvers: optimal maneuver times (s) 
        - manuevers: maneuvers to be executed at t_manuevers (m/s)
        - total_cost: total cost of the reconfiguration (m/s)
    %}
    %% Determine the change in ROE required
    % Make sure you account for any changes to the ROEs that may come from
    % the drift in J2! (STM may be helpful here)
    a = chief_oe(1);
    desired_roe = desired_aroe / a;
    initial_roe = initial_aroe / a;
    Droe = desired_roe - STM(chief_oe, dt) * initial_roe;
    aDroe = a * Droe; 

    %% Calculate the drift rates
    % Here we want to calculate the drift rates of RAAN, AOP, and M
    n   = sqrt(mu/chief_oe(1)^3);                  % mean motion
    eta   = sqrt(1 - chief_oe(2)^2);
    kappa = 3/4 * J2*Rp^2*sqrt(mu) / (chief_oe(1)^(7/2)*eta^4);
    P     = 3*cos(chief_oe(3))^2 - 1;
    Q     = 5*cos(chief_oe(3))^2 - 1;

    RAAN_dot = -2*cos(chief_oe(3))*kappa;
    aop_dot  = kappa * Q;
    M_dot    = n + kappa*eta*P;

    %% Calculate the reachable delta-v
    % Implement calculate_reachable_delta_v()
    [dvmin_da, dvmin_dlambda, dvmin_de, dvmin_di] = calculate_reachable_delta_v(aDroe, chief_oe, dt, n, STM, CM);

    % Then make sure you determine the IP dominance case. To make things 
    % easier, also store the index of the dominance case (look into the max 
    % function from MATLAB, it can return two values!). 
    % 
    % The convention here is: dvmin_dom == 1 - dominance da, 
    %                         dvmin_dom == 2 - dominance dlambda, 
    %                         dvmin_dom == 3 - dominance de
    [dvmin, dvmin_dom] = max([dvmin_da, dvmin_dlambda, dvmin_de]);

    %% Compute the optimal maneuver times
    % Implement compute_optimal_maneuver_times_ip(), and
    % compute_optimal_maneuver_times_oop()

    [T_opts_ip, m_ip] = compute_optimal_maneuver_times_ip(aDroe, chief_oe, dt, n, aop_dot, M_dot);
    [T_opts_oop, m_oop] = compute_optimal_maneuver_times_oop(aDroe, chief_oe, dt, n, kappa, eta, P, Q);

    %% Form the nested reachable sets
    % Following Algorithm 6.1, for each of the optimal times you calculated
    % above, determine the nested reachable set and delta v. 
    
    [Sn_ip, dv_ip] = compute_nested_reachable_sets_and_deltav_ip(T_opts_ip, m_ip, aDroe, chief_oe, dt, dvmin, dvmin_dom, STM, CM, RAAN_dot, aop_dot, M_dot);
    [Sn_oop, dv_oop] = compute_nested_reachable_sets_and_deltav_oop(T_opts_oop, m_oop, chief_oe, dt, dvmin_di, STM, CM, RAAN_dot, aop_dot, M_dot);

    %% Generate all maneuver schemes and pick the best one
    % Implement determine_best_maneuver_plan_ip(), and
    % determine_best_maneuver_plan_oop() using solve_linear_system.
    
    % NOTE: manually edited from aDroe to Droe!
    [ts_ip, dvs_ip, cost_ip] = determine_best_maneuver_plan_ip(Sn_ip, dv_ip, T_opts_ip, m_ip, Droe, dvmin, dvmin_dom);
    % [ts_oop, dvs_oop, cost_oop] = determine_best_maneuver_plan_oop(dv_oop, T_opts_oop, dvmin_di);
  
    %% Create the maneuver plan
    % You're all done! Now we just need to sort the maneuvers so they are
    % returned in the order they occur. Consider the sort function!
    % Don't forget to calculate the total cost!

    % Maneuver and times
    t_maneuvers = ts_ip
    manuevers = dvs_ip

    % Reconfig cost
    total_cost = cost_ip
end

function [dvmin_da, dvmin_dlambda, dvmin_de, dvmin_di] = calculate_reachable_delta_v(aDroe, chief_oe, dt, n, STM, CM)
    %{
    calculate_reachable_delta_v calculates the reachable delta v according
    to Chernick Table 5.13 and 5.14. 

    - Make sure you properly calculate Delta delta lambda_0 and a_0 for the 
        delta lambda dominance case. 
    
    Inputs:
        - aDroe: desired change in dimensionaized QNS ROE state
                a_c Delta [da, dlambda, dex, dey, dix, diy]
        - chief_oe: chief's keplerian orbit elements [a, e, i, RAAN, w, M] 
        - dt: reconfiguration window time (s)
        - n: mean motion of the chief
        - STM: function which given the chief OE's, and a time will return
                the STM by which to multiply ROEs to get the current ROE 
                state
        - CM: returns the control input matrix given the chief's ROEs
        
    Outputs: 
        - dvmin_da: minimum dv for da reconfiguration (m/s)
        - dvmin_dlambda: minimum dv for dlambda reconfiguration (m/s)
        - dvmin_de: minimum dv for de reconfiguration (m/s)
        - dvmin_di: minimum dv for di reconfiguration (m/s)
    %}
    % unpack AOE
    a = chief_oe(1);
    e = chief_oe(2);
    i = chief_oe(3);
    Omega = chief_oe(4);
    w = chief_oe(5);

    % dvmin - dom da
    aDda_des = aDroe(1);
    aDda_max = 2/n;

    dvmin_da = abs(aDda_des)/abs(aDda_max);

    % dmin - dom de
    aDde_max_norm = 2/n;

    aDde_des = aDroe(3:4);
    % page 51 in Chernick (how is this not a dedicated equation!)
    % update: so you're telling me I didn't need to figure this de~ shit
    % out because we take a norm at the end of it all.
    % guh.
    % R = [[cos(w), -sin(w)];
    %      [sin(w), cos(w)]];
    % Dde_tilde = R.' * Dde_des;
    % 
    % aDde_des = a * Dde_tilde;

    dvmin_de = norm(aDde_des) / aDde_max_norm; 

    % dvmin - dom dlambda - limit to tangential burns only
    % no idea if this is tangential-only...
    Dda_des = aDda_des / a;
    aDdlambda_des = aDroe(2);
    Ddlambda_des = aDdlambda_des / a;
    modified_oe = [chief_oe(1:5); 0];
    final_perigee_time = (n*dt - mod(n*dt, 2*pi)) / n;
    
    % STM: t initial is first perigee time, tfinal is last periapsis time
    pseudostate = STM(modified_oe, final_perigee_time) * CM(modified_oe) * [0; 1; 0];
    
    aDda_0 = pseudostate(1);
    aDdlambda_0 = pseudostate(2);

    m = -2 * abs(aDda_0) / abs(aDdlambda_0);

    % TODO: revisit?
    dvmin_dlambda = n*a*abs((m*Ddlambda_des - Dda_des)/aDda_0);

    % dvmin - dom di
    aDdi_max_norm = 1/n;
    aDdi_des = aDroe(5:6);
    % aDdi_tilde = R.' * aDdi_des; % silly me
    dvmin_di = norm(aDdi_des) / aDdi_max_norm;

end

function [T_opts_ip, m_ip] = compute_optimal_maneuver_times_ip(aDroe, chief_oe, dt, n, aop_dot, M_dot)
    %{
    compute_optimal_maneuver_times_ip the optimal maneuver times for in
    plane reconfigurations. Section 6.4 is applied here in the
    near-circular case. Equations 6.9 and 6.10 are particularly helpful. 

    Careful! First define the optimal time with k=m=0. Then, make sure that
    this first maneuver isn't in the past (ie t < 0) by incremeneting m. 
    Once you have a maneuver in the future, figure out at which times 
    within your maneuver window you could complete this maneuver (by 
    increasing k). The end result should be a  list of maneuver locations 
    and times within your maneuver window at which it is optimal to 
    maneuver.
    
    Inputs:
        - aDroe: desired change in dimensionaized QNS ROE state
                a_c Delta [da, dlambda, dex, dey, dix, diy]
        - chief_oe: chief's keplerian orbit elements [a, e, i, RAAN, w, M] 
        - dt: reconfiguration window time (s)
        - n: mean motion of the chief
        - aop_dot: time rate of change of the argument of periapsis under
                J2
        - M_dot: time rate of change of the mean anomaly under J2
        
    Outputs: 
        - T_opts_ip: all the maneuver times at which we could optimally
                maneuver to affect the desired in plane ROEs
        - m_ip: value of m needed to get the first maneuver time to be in
                the future (ie t > 0)
    %}
    % TODO Calculate the optimal time to maneuver (may be in the past!)
    a = chief_oe(1);
    w = chief_oe(5);

    Ddex = aDroe(3) / a;
    Ddey = aDroe(4) / a;
    
    % compute m:
    if atan2(Ddey,Ddex) - (aop_dot*dt + w) < 0
        m_ip = ceil((aop_dot*dt + w) - atan2(Ddey,Ddex) / pi);
    else
        m_ip = 0;
    end

    % compute k
    T0 = (atan2(Ddey,Ddex) + m_ip*pi - (aop_dot*dt + w))/M_dot;
    k_max = floor((dt - T0) * M_dot / pi);

    T_opts_ip = T0 + (0:k_max) * pi / M_dot;
    
end

function [T_opts_oop, m_oop] = compute_optimal_maneuver_times_oop(aDroe, chief_oe, dt, n, kappa, eta, P, Q)
    %{
    compute_optimal_maneuver_times_oop the optimal maneuver times for out
    of plane reconfigurations. Section 6.4 is applied here in the
    near-circular case. Equation 6.11 is particularly helpful. 

    Careful! First define the optimal time with k=m=0. Then, make sure that
    this first maneuver isn't in the past (ie t < 0) by incrementing m. 
    Once you have a maneuver in the future, figure out at which times 
    within your maneuver window you could complete this maneuver (by 
    increasing k). The end result should be a  list of maneuver locations 
    and times within your maneuver window at which it is optimal to 
    maneuver.
    
    Inputs:
        - aDroe: desired change in dimensionaized QNS ROE state
                a_c Delta [da, dlambda, dex, dey, dix, diy]
        - chief_oe: chief's keplerian orbit elements [a, e, i, RAAN, w, M] 
        - dt: reconfiguration window time (s)
        - n: mean motion of the chief
        - kappa, eta, P, Q: constants defined in eqn 2.6
        
    Outputs: 
        - T_opts_oop: all the maneuver times at which we could optimally
                maneuver to affect the desired out of plane ROEs
        - m_oop: value of m needed to get the first maneuver time to be in
                the future (ie t > 0)
    %}
    % Calculate the optimal time to maneuver (may be in the past!)
    a = chief_oe(1);
    w = chief_oe(5);

    Ddix = aDroe(5) / a;
    Ddiy = aDroe(6) / a;
    
    denom = n + kappa*(eta*P+Q);

    % Add to m until the time at this desired phase is greater than 0
    m_oop = ceil((1/pi)*(w-atan2(Ddiy, Ddix)));
    numerator = atan2(Ddiy, Ddix) - w + m_oop*pi;
    T0 = numerator / denom;

    % Determine when are all the times we could do this maneuver within our
    % maneuver window
    k_max = floor((denom*(dt - T0) - numerator)/pi);

    T_opts_oop = T0 + (0:k_max) * pi/denom;
end

function [Sn_ip, dv_ip] = compute_nested_reachable_sets_and_deltav_ip(T_opts_ip, m_ip, aDroe, chief_oe, dt, dvmin, dvmin_dom, STM, CM, RAAN_dot, aop_dot, M_dot)
    %{
    compute_nested_reachable_sets_and_deltav_ip calculates the nested
    reachable sets for the in plane reconfiguration and the delta v unit
    vectors. These delta v's should be unitless, no need to dimensionalize
    them yet. 

    - To determine the delta v, refer to Table 6.3 
    - The definition of the nested reachable set is on Algorithm 6.1 
        (line 6). 
    - When creating the nested reachable set, make sure you use the chief 
        orbit elements at the right time. Can you think of a simple way to 
        estimate the chief orbit elements at any time under J2 perturbation? 
    - For the dlambda dominance, see table 7.1.
    
    Inputs:
        - T_opts_ip: all the maneuver times at which we could optimally
                maneuver to affect the desired in plane ROEs
        - m_ip: value of m needed to get the first maneuver time to be in
                the future (ie t > 0)
        - aDroe: desired change in dimensionaized QNS ROE state
                a_c Delta [da, dlambda, dex, dey, dix, diy]
        - chief_oe: chief's keplerian orbit elements [a, e, i, RAAN, w, M] 
        - dt: reconfiguration window time (s)
        - dvmin: delta v minimum of the dominance case (m/s)
        - dvmin_dom: dominance case (1: da, 2: dlambda, 3: de)
        - STM: function which given the chief OE's, and a time will return
                the STM by which to multiply ROEs to get the current ROE 
                state
        - CM: returns the control input matrix given the chief's ROEs
        - RAAN_dot: time rate of change of the right ascencion of the
                ascending node under J2 
        - aop_dot: time rate of change of the argument of periapsis under
                J2
        - M_dot: time rate of change of the mean anomaly under J2
        
    Outputs: 
        - Sn_ip - nested reachable sets
        - dv_ip - delta v's UNIT VECTORS
    %}
    Sn_ip = zeros(6, length(T_opts_ip)); % nested reachable sets
    dv_ip = zeros(3, length(T_opts_ip));

    for k = 1:length(T_opts_ip)
        %  delta v
        if dvmin_dom == 1 % dom da
            if aDroe(1) > 1
                dv_ip(:, k) = [0; 1; 0];
            else
                dv_ip(:, k) = [0; -1; 0];
            end
        elseif dvmin_dom == 2 % dom dlambda
            % what's going on in Table 7.1??
            if aDroe(2) < 0 && k == 1 || aDroe(2) > 0 && k==length(T_opts_ip)
                dv = 1;
            elseif (aDroe(2) > 0 && k == 1) || (aDroe(2) < 0 && k==length(T_opts_ip))
                dv = -1;
            else
                if mod(m_ip+k-3, 2) == 0
                    dv = 1;
                else
                    dv = -1;
                end
            end
            dv_ip(:, k) = [0; dv; 0];
        elseif dvmin_dom == 3 % dom de
            if mod(k+m_ip, 2) == 0
                dv_ip(:, k) = [0; 1; 0];
            else
                dv_ip(:, k) = [0; -1; 0];
            end
        end

        % Nested reachable set
        % calculate chief OE at maneuver time
        t = T_opts_ip(k);
        % Runge and Kutta are rolling in their graves
        J2_propagation = wrapTo2Pi([0, 0, 0, RAAN_dot*t, aop_dot*t, M_dot*t]);
        % [a, e, i, RAAN, w, M]
        chief_oe_k = chief_oe + J2_propagation;
        
        Gamma = STM(chief_oe_k, dt - t) * CM(chief_oe_k);
        Sn_ip(:,k) = Gamma * dv_ip(:,k) * dvmin;
    end
end

function [Sn_oop, dv_oop] = compute_nested_reachable_sets_and_deltav_oop(T_opts_oop, m_oop, chief_oe, dt, dvmin_di, STM, CM, RAAN_dot, aop_dot, M_dot)
    %{
    compute_nested_reachable_sets_and_deltav_oop calculates the nested
    reachable sets for the out of plane reconfiguration and the delta v 
    unit vectors. These delta v's should be unitless, no need to dimensionalize
    them yet. 

    - To determine the delta v, refer to Table 6.3 
    - The definition of the nested reachable set is on Algorithm 6.1 
        (line 6). 
    - When creating the nested reachable set, make sure you use the chief 
        orbit elements at the right time. Can you think of a simple way to 
        estimate the chief orbit elements at any time under J2 perturbation? 
    
    Inputs:
        - T_opts_oop: all the maneuver times at which we could optimally
                maneuver to affect the desired out of plane ROEs
        - m_oop: value of m needed to get the first maneuver time to be in
                the future (ie t > 0)
        - aDroe: desired change in dimensionaized QNS ROE state
                a_c Delta [da, dlambda, dex, dey, dix, diy]
        - dt: reconfiguration window time (s)
        - dvmin_di: delta v minimum of the dominance case (di always) (m/s)
        - STM: function which given the chief OE's, and a time will return
                the STM by which to multiply ROEs to get the current ROE 
                state
        - CM: returns the control input matrix given the chief's ROEs
        - RAAN_dot: time rate of change of the right ascencion of the
                ascending node under J2 
        - aop_dot: time rate of change of the argument of periapsis under
                J2
        - M_dot: time rate of change of the mean anomaly under J2
        
    Outputs: 
        - Sn_ip - nested reachable sets
        - dv_ip - delta v's
    %}
    Sn_oop = zeros(6, length(T_opts_oop)); % nested reachable sets
    dv_oop = zeros(3, length(T_opts_oop));

    for k = 1:length(T_opts_oop)
        % dom di - always
        km = k + m_oop;
        if mod(km, 2) == 0
            dv_oop(:,k) = [0; 0; 1];
        else
            dv_oop(:,k) = [0; 0; -1];
        end
        

        % Nested reachable set
        % calculate chief OE at maneuver time
        t = T_opts_oop(k);
        % Runge and Kutta are rolling in their graves
        J2_propagation = wrapTo2Pi([0, 0, 0, RAAN_dot*t, aop_dot*t, M_dot*t]);
        % [a, e, i, RAAN, w, M]
        chief_oe_k = chief_oe + J2_propagation;
        
        Gamma = STM(chief_oe_k, dt - t) * CM(chief_oe_k);
        Sn_oop(:, k) = Gamma * dv_oop(:, k) * dvmin_di;
    end
end

function [cs] = solve_linear_system(M,b)
    if rcond(M) < 1e-6
        cs = [-1;-1;-1];
    else
        cs = M \ b;
    end
end

function [ts_ip, dvs_ip, cost_ip] = determine_best_maneuver_plan_ip(Sn_ip, dv_ip, T_opts_ip, m_ip, aDroe, dvmin, dvmin_dom)
    %{
    determine_best_maneuver_plan_ip chooses the best maneuver plan for the
    in plane reconfiguration which minimizes the cost.

    Please see Table 6.4. This is likely the trickiest function to
    implement, so be careful. We want to set up a linear system which is
    defined row-wise on the third column of this table. 
    
    Make sure that in the da dominant case, you take into account the 
    effect of maneuver on de vector wrt desired pseudostate.

    Note that there is also the sub-optimal case, where instead of the top
    row being 1's, it is the de.

    We include the function `solve_linear_system` to help you handle
    numerical edge cases. 
    
    Inputs:
        - Sn_ip - nested reachable sets
        - dv_ip - delta v's
        - T_opts_ip: all the maneuver times at which we could optimally
                maneuver to affect the desired in plane ROEs
        - m_ip: value of m needed to get the first maneuver time to be in
                the future (ie t > 0)
        - aDroe: desired change in dimensionaized QNS ROE state
                a_c Delta [da, dlambda, dex, dey, dix, diy]
        - dvmin: delta v minimum of the dominance case (m/s)
        - dvmin_dom: dominance case (1: da, 2: dlambda, 3: de)
        
    Outputs: 
        - ts_ip - optimal maneuver times (s)
        - dvs_ip - maneuvers to be exectuted at the optimal maneuver times
            now dimensionalized (m/s)
        - cost_ip - total cost of the reconfiguration (m/s)
    %}
    
    dvs_ip = [];
    ts_ip = [];
    cost_ip = inf;
    combos = nchoosek(1:length(T_opts_ip), 3);

    for i = 1:size(combos,1)
        inds = combos(i,:);

        % Nested reachable sets
        viable_Sn = Sn_ip(:, inds);
        
        if dvmin_dom == 1 % dom da
            % TODO Effect of maneuver on de vector wrt desired pseudostate

            % TODO Optimal linear system
            cs = 0.5;
            
        elseif dvmin_dom == 2 % dom dl
            % TODO Optimal linear system - none
            % cs = ??;
            cs = 0.5;
        elseif dvmin_dom == 3 % dom de
            % TODO Optimal linear system - Table 6.4

            Dda1 = viable_Sn(:, 1);
            Dda2 = viable_Sn(:, 2);
            Dda3 = viable_Sn(:, 3);

            A = [1,        1,         1;
                 Dda1(1),  Dda2(1),   Dda3(1);
                 Dda1(2),  Dda2(2),   Dda3(2);];
            b = [1; aDroe(1); aDroe(2)];
            cs = solve_linear_system(A, b);
        end 

        if any(cs < 0) || any(cs > 1) % sub optimal solution
            % TODO Effect of maneuver on de vector wrt desired pseudostate
            
            % TODO Suboptimal linear system
            A = [norm(Dda1(3:4)),  norm(Dda2(3:4)),   norm(Dda3(3:4));
                 Dda1(1),          Dda2(1),           Dda3(1);
                 Dda1(2),          Dda2(2),           Dda3(2);];
            b = [norm(aDroe(3:4)); aDroe(1); aDroe(2)];
            cs = solve_linear_system(A, b);
        end

        % TODO Determine the best solution
        dvs_ip_cand = cs .* dv_ip(:, inds);

        if sum(vecnorm(dvs_ip_cand)) < cost_ip
            disp(inds)
            dvs_ip = cs .* dv_ip(:, inds);
            ts_ip = T_opts_ip(inds);
            cost_ip = sum(vecnorm(dvs_ip));
        end
    end

end

function [ts_oop, dvs_oop, cost_oop] = determine_best_maneuver_plan_oop(dv_oop, T_opts_oop, dvmin_di)
    %{
    determine_best_maneuver_plan_oop chooses the best maneuver plan for the
    out of plane reconfiguration which minimizes the cost.

    Note that in the out of plane cases, all maneuvers have the same cost.
    Feel free to just pick the last maneuver that is possible in the
    reconfiguration window. 
    
    Inputs:
        - dv_oop - delta v's
        - T_opts_oop: all the maneuver times at which we could optimally
                maneuver to affect the desired out of plane ROEs
        - - dvmin_di: delta v minimum of the dominance case (di always) (m/s)
        
    Outputs: 
        - ts_ip - optimal maneuver times (s)
        - dvs_ip - maneuvers to be exectuted at the optimal maneuver times
            now dimensionalized (m/s)
        - cost_ip - total cost of the reconfiguration (m/s)
    %}
    
    % TODO
    % dvs_oop = ??;
    % ts_oop = ??;
    % cost_oop = ??;
end
