clc
clear all
close all
%% Technical parameters
options    = optimoptions('fsolve', 'Display', 'iter', 'OptimalityTolerance', 1e-6);
options_c  = optimoptions('fsolve', 'Display', 'none', 'OptimalityTolerance', 1e-6);
%% Parameters
beta = 0.97;              % Discount factor
alpha = 0.33;             % Capital share
psi   = 1  ;              % Frisch Elasticity
delta = 0.1;              % Depreciation rate
rho = 0.95;               % AR(1) persistence of z_t
sigma_eps = 0.007;        % Std. deviation of shocks
num_z = 3;                % Number of grid points for z_t (Tauchen's method)
num_cheb = 6;             % Number of Chebyshev polynomials
z_ss = 0;
%% Step 1: Discretize z_t using Tauchen's method
[z_grid,P] = MC_Tauchen(num_z,z_ss,rho,sigma_eps,3);
z_grid = exp(z_grid');

%% STEP 1.1 SP Steady State             
SP_steady_syst = @(x) [(1/beta) - (alpha*x(2)^(alpha-1)*x(3)^(1-alpha)+1-delta);
                        delta*x(2) - x(2)^(alpha)*x(3)^(1-alpha) + x(1); 
                        x(3)^psi - (1-alpha)*x(2)^(alpha)*x(3)^(-alpha)*(1/x(1));
                    ];
                
        
x0 = ones(3,1);
[xce_sol,fval,exitflag,output] = fsolve(SP_steady_syst,x0,options);                

c_ss = xce_sol(1);
k_ss = xce_sol(2);
l_ss = xce_sol(3);


%% Step 2: Define Chebyshev nodes for capital
k_min = k_ss*0.75;              % Minimum capital
k_max = k_ss*1.25;                % Maximum capital
k_grid = linspace(k_min, k_max, num_cheb);
k_nodes = cos((2 * (1:num_cheb) - 1) / (2 * num_cheb) * pi);  % Chebyshev zeros
% k_nodes = k_min + (k_max - k_min) * (cheb_nodes + 1) / 2;        % Transform to [k_min, k_max]
% k_nodes = k_grid;

%% Step 3: Initialize coefficients for Chebyshev approximation
theta_k = ones(num_cheb, num_z)/(num_z*num_cheb);  % Coefficients for each z_t level

%% Step 4: Define labor and consumption functions
% Labor as a function of consumption
labor = @(c, k, z) ((1 - alpha)*z.* k.^alpha ./ c).^(1 / (psi + alpha));

% Solve for consumption using the resource constraint
solve_consumption = @(c,k,kp,z) z*k^alpha*labor(c, k, z)^(1 - alpha)-c - kp +(1-delta)*k;

%% Step 5: Residual function for Euler equation
residual = @(theta_k) deal(residual_function(theta_k, k_nodes,k_grid,z_grid, P, beta, alpha, delta, labor, solve_consumption,options_c));

%% Step 6: Solve for coefficients
tic
theta_k = fsolve(residual, theta_k(:), options);
toc
%% Step 7: Post-process and plot results
theta_k = reshape(theta_k, [num_cheb, num_z]);  % Reshape coefficients
plot_policy_functions(k_nodes, k_grid,z_grid, theta_k,solve_consumption,options_c);
