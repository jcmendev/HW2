function res = residual_function_fem(theta_k, k_grid, z_grid,k_nodes, P, beta, alpha, delta, labor, solve_consumption, basis_function,options_c)
    num_k = length(k_grid);
    num_z = length(z_grid);
    res = zeros(num_k, num_z);

    mtheta_k = reshape(real(theta_k),num_k,num_z);
   
    for iz = 1:num_z
        z = z_grid(iz);

        for ik = 1:num_k
            k = k_grid(ik);
            % Approximate next-period capital using finite elements
            kp_approx = 0;
            for i = 1:num_k
                kp_approx = kp_approx + mtheta_k(i, iz) * basis_function(k, k_grid(max(i-1, 1)), k_grid(min(i+1, num_k)));
            end
            
%             kpm = kp_approx;
            km  = (((k+1)/2)*(k_nodes(end)-k_nodes(1)))+k_nodes(1);
            kpm = (((kp_approx+1)/2)*(k_nodes(end)-k_nodes(1)))+k_nodes(1);
            
%             km1 = (k*(k_nodes(min(ik+1, num_k))-k_nodes(max(ik-1, 1)))+(k_nodes(max(ik-1, 1))+k_nodes(min(ik+1, num_k))))/2;
%             kpm1 = (kp_approx*(k_nodes(min(ik+1, num_k))-k_nodes(max(ik-1, 1)))+(k_nodes(max(ik-1, 1))+k_nodes(min(ik+1, num_k))))/2;

%             kpm = ((kp_approx*(k_grid(ik+1)-k_grid(ik-1)))+k_grid(ik-1);            
%             kp_approx = min(max(kp_approx,k_grid(1)),k_grid(end));
            % Solve for current-period consumption
            
            residual_c = @(c) deal(solve_consumption(c,km,kpm,z));
            c          = (fsolve(residual_c, 1, options_c));            

            % Labor decision
            l = labor(c, k, z);

            % Compute Euler equation residual
            expected_val = 0;
            for jz = 1:num_z
                z_next = z_grid(jz);
                kp_next = 0;
                
                kp_approx_aux = 2*(kpm-k_nodes(1))/(k_nodes(end)-k_nodes(1))-1 ;
                for i = 1:num_k
                    kp_next = kp_next + mtheta_k(i, jz) * basis_function(kp_approx, k_grid(max(i-1, 1)), k_grid(min(i+1, num_k)));
                end
                
%                 kp_next_m = (((kp_next)/2)*(k_grid(end)-k_grid(1)))+k_grid(1);
                kp_next_m = (((kp_next+1)/2)*(k_nodes(end)-k_nodes(1)))+k_nodes(1);
%                 kp_next_m = (kp_next*(k_nodes(min(ik+1, num_k))-k_nodes(max(ik-1, 1)))+(k_nodes(max(ik-1, 1))+k_nodes(min(ik+1, num_k))))/2;
                
%                 kp_next = min(max(kp_next,k_grid(1)),k_grid(end));
                residual_c_next = @(cp) deal(solve_consumption(cp,kpm,kp_next_m,z_next));
                c_next          = (fsolve(residual_c_next, 1, options_c));                
                l_next = labor(c_next, kpm, z_next);

                mpk = alpha*z_next*kp_next_m^(alpha - 1)*l_next^(1 - alpha) + (1 - delta);
                expected_val = expected_val + P(iz, jz) * (1 / c_next) * mpk;
            end

            % Galerkin weighting: integrate residual against basis functions
            res_aux = 0;
            for ikk = 1:num_k
                res_aux = res_aux + (basis_function(k, k_grid(max(ikk-1, 1)), k_grid(min(ikk+1, num_k))));
            end
            res(ik, iz) = res_aux*(1/c - beta*expected_val);
        end
    end

    res = res(:);  % Flatten residuals into a vector
end
