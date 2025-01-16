%% Variables
var
        c
        k
        l
        z         
;

varexo 
sh_z_t
;

%% Parameters
parameters

beta 
alpha 
psi   
delta 
rho 
sigma_eps 
;

beta = 0.97;              % Discount factor
alpha = 0.33;             % Capital share
psi   = 1  ;              % Frisch Elasticity
delta = 0.1;              % Depreciation rate
rho = 0.95;               % AR(1) persistence of z_t
sigma_eps = 0.007; 

%% MODEL
model; 
%% System of equations 
% % ---- Eq 1: NAME ---- % %
exp(l)^psi = (1-alpha)*exp(z)*exp(k(-1))^alpha*exp(l)^(-alpha)*(1/exp(c));
% % ---- Eq 2: NAME ---- % %
(1/exp(c)) = beta*(((alpha)*exp(z(+1))*exp(k)^(alpha-1)*exp(l(+1))^(1-alpha)+1-delta))*(1/exp(c(+1)));
% % ---- Eq 3: NAME ---- % %
exp(k) = (1-delta)*exp(k(-1)) + exp(z)*exp(k(-1))^alpha*exp(l)^(1-alpha) - exp(c);
% % ---- Eq 4: NAME ---- % %
z   =  rho*(z(-1)) + sh_z_t;

end;

%% Steady state
initval; 
z = 0;
c = log(1.1161);
l = log(0.9465);
k = log(3.7612);
end;
steady;
resid;

%% Third Order Perturbation
shocks;
var sh_z_t;
stderr sigma_eps;
end;
stoch_simul (order=3, pruning,irf=0);

