function plot_policy_functions(k_nodes, k_grid,z_grid, theta_k,solve_consumption,options_c)
    num_k = length(k_nodes);
    num_z = length(z_grid);
    num_cheb = length(k_nodes);
    for iz = 1:num_z
        z = z_grid(iz);
        for ik = 1:num_cheb
            k = k_nodes(ik);
            kp(ik,iz) = sum(theta_k(:, iz) .* chebyshev_polynomials(k, k_nodes));
        end
    end
    
    km  = (((k_nodes+1)/2)*(k_grid(end)-k_grid(1)))+k_grid(1);
    kpm = (((kp+1)/2)*(k_grid(end)-k_grid(1)))+k_grid(1);
            
    figure;
    for iz = 1:num_z
        z = z_grid(iz);
        plot(km, kpm(:,iz), 'DisplayName', ['z = ', num2str(z)]);
        hold on;
    end
    plot(km,km,'--k')
    xlabel('Capital (k)','Interpreter','Latex');
    ylabel('Next-period Capital (k'')','Interpreter','Latex');
    title('Policy Functions using Chebyshev Approximation','Interpreter','Latex');
    legend('show');
    grid on;
    print(['Fig_Cheb_K'],'-depsc','-r0')


% WARNING OFF
    for ik = 1:num_k
    for iz = 1:num_z
        z = z_grid(iz);
        residual_c = @(c) deal(solve_consumption(c,km(ik),kpm(ik,iz),z));
        c(ik,iz)   = (fsolve(residual_c, 1, options_c));

    end
    end
    
        figure;
    for iz = 1:num_z    
        z = z_grid(iz);
        plot(km, c(:,iz), 'DisplayName', ['z = ', num2str(z)]);
        hold on;
    end
    xlabel('Capital (k)','Interpreter','Latex');
    ylabel('Consumption','Interpreter','Latex');
    title('Policy Functions using Chebyshev Approximation','Interpreter','Latex');
    legend('show');
    grid on;
    print(['Fig_Cheb_C'],'-depsc','-r0')    

end
