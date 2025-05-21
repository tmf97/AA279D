function continuous_control_wrapper(roes_desired, mu, N, k, phi_ip, phi_oop)
    function ds = continuous_control(~, s)
        % compute Keplerian and J2 accelerations for both chief and deputy
        r_chief = s(1:3);
        v_chief = s(4:6);
        r_deputy = s(7:10);
        v_deputy = s(11:13);
        
        a_chief = newton3d_J2(0, [r_chief; v_chief]);
        a_deputy = newton3d_J2(0, [r_deputy; v_deputy]);
        
        % compute KOEs for both chief and deputy
        koe_chief = pv2koe([r_chief; v_chief], mu);
        koe_deputy = pv2koe([r_deputy; v_deputy], mu);
        
        moe_chief = osc2mean(koe_chief, 1);
        moe_deputy = osc2mean(koe_deputy, 1);
        roes = quasi_nonsingular_roe(moe_chief, moe_deputy);
        ddalpha = roes - roes_desired;
        
        e = a_deputy(2);
        nu = a_deputy(6);
        w = a_deputy(5);
        M = true2mean(nu, e);

        % mean argument of latitude of the deputy
        phi = w + M;

        cosJ = cos(phi - phi_ip)^N;
        cosH = cos(phi - phi_oop)^N;
        
        P = (1/k) * diag([cosJ, cosJ, cosJ, cosJ, cosH, cosH]);
        A = plant_matrix_J2(koe_chief, mu);
        B = control_input_matrix_v2(koe_chief, mu);
        Binv = pinv(B);

        u_RTN = -Binv * (A * roes + P * ddalpha);
        
        rhat = r_last/norm(r_chief);
        h0 = cross(r_chief, v_chief);
        nhat = h0/norm(h0);
        that = cross(nhat, rhat);
        C_RTN_to_ECI = [rhat, that, nhat];

        u_ECI = C_RTN_to_ECI * u_RTN;
        
        ds = [v_chief; v_deputy; a_chief; a_deputy+u_ECI];
    end
end