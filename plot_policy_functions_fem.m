function plot_policy_functions_fem(k_grid, z_grid, theta_k, basis_function,k_nodes)

    num_z = length(z_grid);
    num_k = length(k_grid);
    kp = zeros(num_k,num_z);
    for iz = 1:num_z
        z = z_grid(iz);
        
        for ik = 1:num_k
            k = k_grid(ik);
            kp(ik,iz) = 0;
            for i = 1:num_k
                kp(ik,iz) = kp(ik,iz) + theta_k(i, iz) * basis_function(k, k_grid(max(i-1, 1)), k_grid(min(i+1, length(k_grid))));
            end
        end

    end
    
    for iz = 1:num_z    
    for ik = 1:num_k
        k = k_grid(ik);
        kp_approx = kp(ik,iz);
            km(ik,iz) = (k*(k_nodes(min(ik+1, num_k))-k_nodes(max(ik-1, 1)))+(k_nodes(max(ik-1, 1))+k_nodes(min(ik+1, num_k))))/2;
            kpm(ik,iz) = (kp_approx*(k_nodes(min(ik+1, num_k))-k_nodes(max(ik-1, 1)))+(k_nodes(max(ik-1, 1))+k_nodes(min(ik+1, num_k))))/2;
    end
    end            
    km  = (((k_grid+1)/2)*(k_nodes(end)-k_nodes(1)))+k_nodes(1);
    kpm = (((kp+1)/2)*(k_nodes(end)-k_nodes(1)))+k_nodes(1);
%     kpm = kp;
    figure;  
    for iz = 1:num_z   
        z = z_grid(iz);
        plot(km, kpm(:,iz), 'DisplayName', ['z = ', num2str(z)]);
        hold on;
    end     
    plot(km,km,'--k')
    xlabel('Capital (k)','Interpreter','Latex');
    ylabel('Next-period Capital (k'')','Interpreter','Latex');
    title('Policy Functions using Finite Elements Approximation','Interpreter','Latex');
    legend('show');
    grid on;
    print(['Fig_FEM_K'],'-depsc','-r0')
end
