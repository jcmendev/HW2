clc
clear all
% close all
%% Addpath
% addpath C:\dynare\4.4.3\matlab
%% Step 1: Discretize z_t using Tauchen's method
z_ss = 0;
rho = 0.95;               % AR(1) persistence of z_t
sigma_eps = 0.007;        % Std. deviation of shocks
num_z = 3;                % Number of grid points for z_t (Tauchen's method)

[z_grid,P] = MC_Tauchen(num_z,z_ss,rho,sigma_eps,3);
z_grid = z_grid';
% z_grid = exp(z_grid);

%% Step 2: Run Dynare
tic
dynare('Perturbation_Dynare.mod','noclearall')
toc
%% Step 3: Load Dynare output

load Perturbation_Dynare_results.mat
ghx = oo_.dr.ghx;  % Coefficients on state variables
ghu = oo_.dr.ghu;  % Coefficients on shocks
steady_state = oo_.steady_state;  % Steady-state values of variables

c_ss = steady_state(find(strcmp(cellstr(M_.endo_names), 'c')));  % Steady-state capital
k_ss = steady_state(find(strcmp(cellstr(M_.endo_names), 'k')));  % Steady-state capital
l_ss = steady_state(find(strcmp(cellstr(M_.endo_names), 'l')));  % Steady-state capital

%% Step 2: Define grid for capital
num_k = 100;                % Number of finite elements for capital
k_min = k_ss*0.75;              % Minimum capital
k_max = k_ss*1.25;                % Maximum capital
k_grid    = linspace(k_min, k_max, num_k);  % Capital grid

% Compute deviations of k_t from steady state
k_dev = k_grid - k_ss;
z_dev = z_grid - z_ss;

%% Extract info
vtp = {'Capital',1,k_ss};

vss = vtp{3};
% First-order coefficients
ghx_kp = oo_.dr.ghx(vtp{2},:);

% Second-order coefficients
ghxx_kp = oo_.dr.ghxx(vtp{2},:);
rghxx_kp = reshape(ghxx_kp,2,2);

% Third-order coefficients
ghxxx_kp = oo_.dr.ghxxx(vtp{2},:);
rghxxx_kp = reshape(ghxxx_kp,2,2,2);
%%

for i = 1:num_k
    for j = 1:num_z
        % Deviations from steady-state
        state_dev = [k_dev(i); z_dev(j)];
        
        % First-order term
        First = vss + ghx_kp * state_dev;
        
        Second = 0.5 * state_dev' * rghxx_kp * state_dev;
        
        Third = 0;
        for in = 1:2
            Third = Third + (1/6) * state_dev'*rghxxx_kp(:,:,in)* state_dev;
        end
        
        policy(i, j) = First  + Second + Third;
    end
end

%%
figure()
    for iz = 1:num_z
        z = z_grid(iz);
        plot(exp(k_grid),exp(policy(:,iz)), 'DisplayName', ['z = ', num2str(z)]);
        hold on
    end
    plot(exp(k_grid),exp(k_grid),'--k')
xlabel('Current Capital \( k_t \)', 'Interpreter', 'latex');
ylabel('Next Period Capital \( k_{t+1} \)', 'Interpreter', 'latex');
title('Policy Function (Third-Order Approximation)', 'Interpreter', 'latex');
legend('show');
grid on;
print(['Fig_Pert_K'],'-depsc','-r0')