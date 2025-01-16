function res = residual_function(theta_k, k_nodes,k_grid, z_grid, P, beta, alpha, delta, labor, solve_consumption,options_c)
    num_cheb = length(k_nodes);
    num_z = length(z_grid);
    res = zeros(num_cheb, num_z);

    mtheta_k = reshape(theta_k,num_cheb,num_z);
    for iz = 1:num_z
        z = z_grid(iz);

        for ik = 1:num_cheb
            k = k_nodes(ik);
            kp_approx = sum(mtheta_k(:, iz) .* chebyshev_polynomials(k, k_nodes));
            
            km  = (((k+1)/2)*(k_grid(end)-k_grid(1)))+k_grid(1);
            kpm = (((kp_approx+1)/2)*(k_grid(end)-k_grid(1)))+k_grid(1);
            % Solve for current-period consumption
            residual_c = @(c) deal(solve_consumption(c,km,kpm,z));
            c_approx   = (fsolve(residual_c, 1, options_c));

            % Labor decision
            l = labor(c_approx, km, z);

            % Compute Euler equation residual
            expected_val = 0;
            for jz = 1:num_z
                z_next = z_grid(jz);
                kp_next = sum(mtheta_k(:, jz) .* chebyshev_polynomials(kp_approx, k_nodes));
                                
                kp_next_m = (((kp_next+1)/2)*(k_grid(end)-k_grid(1)))+k_grid(1);
            
                residual_cp     = @(cp) deal(solve_consumption(cp,kpm,kp_next_m,z_next));
                c_next   = fsolve(residual_cp, 1, options_c);
                l_next = labor(c_next, kpm, z_next);

                mpk = alpha*z_next*kp_next_m^(alpha - 1)*l_next^(1 - alpha) + (1 - delta);
                expected_val = expected_val + P(iz, jz) * (1 / c_next) * mpk;
            end

            res(ik, iz) = 1/c_approx - beta*expected_val;
        end
    end

    res = res(:);  % Flatten residuals into a vector
end
