clc
clear all
% close all
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
num_k = 6;                % Number of finite elements for capital
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

%% Step 2: Define grid for capital
k_min = k_ss*0.75;              % Minimum capital
k_max = k_ss*1.25;                % Maximum capital
k_grid    = linspace(-1, 1, num_k);  % Capital grid
k_nodes   = linspace(k_min, k_max, num_k);  % Capital grid
% k_grid = k_nodes;
%% Step 3: Initialize coefficients for finite elements
% theta_k = ones(num_k, num_z)/(num_k*num_z);  % Coefficients for each z_t level
theta_k = ones(num_k, num_z)/(num_k*num_z);  % Coefficients for each z_t level

%% Step 4: Define basis functions for finite elements
% Piecewise linear basis functions
basis_function = @(x, xi1, xi2) max(0, (x - xi1) / (xi2 - xi1) .* (x >= xi1 & x <= xi2)) + ...
                                max(0, (xi2 - x) / (xi2 - xi1) .* (x >= xi1 & x <= xi2));

%% Step 5: Define labor and consumption functions
% Labor as a function of consumption
labor = @(c, k, z) ((1 - alpha)*z.* k.^alpha ./ c).^(1 / (psi + alpha));

% Solve for consumption using the resource constraint
solve_consumption = @(c,k,kp,z) z*k^alpha*labor(c, k, z)^(1 - alpha)-c - kp +(1-delta)*k;

%% Step 6: Residual function for Euler equation with Galerkin weighting
residual = @(theta_k) deal(residual_function_fem(theta_k, k_grid, z_grid,k_nodes, P, beta, alpha, delta, labor, solve_consumption, basis_function,options_c));

%% Step 7: Solve for coefficients
tic
theta_k = fsolve(residual, theta_k(:), options);
toc
%% Step 8: Post-process and plot results
theta_k = reshape(theta_k, [num_k, num_z]);  % Reshape coefficients
plot_policy_functions_fem(k_grid, z_grid, theta_k, basis_function,k_nodes);

