clear
clc
close all
format long g
% This is the folder where the VFI toolkit files are saved
folder1 = 'C:\Users\aledi\Desktop\VFIToolkit-matlab-master';
%folder2 = fullfile('..','tools');
addpath(genpath(folder1))

%% Set flags

do_pijoan = 1; % If 1, load shocks from Pijoan-Mas files, otherwise discretize

%% Set economic parameters

par.crra   = 1.458; % Coeff of risk aversion
par.nu     = 2.833; % Curvature labor utility
par.lambda = 0.856; % Weight of labor in util
par.beta   = 0.945; % Discount factor
par.theta  = 0.64;  % Labor share in Cobb-Douglas
par.delta  = 0.083; % Capital depreciation rate

% HSV taxation. 
% No taxes: lambda=1, tau=0
% proportional taxes: lambda=1-taxrate, tau=0
par.lambda_hsv = 1.0;
par.tau_hsv    = 0.0;

% Parameters for AR1 labor productivity z
rho_z = 0.92;
sig_z = 0.21;

%% Set grids and shocks

% Grid for assets
a_curve = 2.0;
a_min   = 1e-6;
a_max   = 50;
n_a     = 601;
a_grid = make_grid(a_min,a_max,n_a,a_curve,1);

%grid for labor
n_d    = 51;
d_grid = linspace(0.001,0.999,n_d)';

if do_pijoan==1

    z_grid_log = [-1.385493 -0.923662 -0.461831  0.000000  0.461831  0.923662  1.385493]';
    z_grid     = exp(z_grid_log);
    pi_z   =  [0.746464  0.252884  0.000652  0.000000  0.000000  0.000000  0.000000
        0.046088  0.761085  0.192512  0.000314  0.000000  0.000000  0.000000
        0.000028  0.069422  0.788612  0.141793  0.000145  0.000000  0.000000
        0.000000  0.000065  0.100953  0.797965  0.100953  0.000065  0.000000
        0.000000  0.000000  0.000145  0.141793  0.788612  0.069422  0.000028
        0.000000  0.000000  0.000000  0.000314  0.192512  0.761085  0.046088
        0.000000  0.000000  0.000000  0.000000  0.000652  0.252884  0.746464];

    pi_z = pi_z./sum(pi_z,2);
    n_z = length(z_grid);

else
    n_z = 7;
    [pi_z,z_grid_log] = markovapprox(rho_z,sig_z,0.0,3.0,n_z,0);
    z_grid     = exp(z_grid_log);

end

% Pack grids into a struct
grids = pack_into_struct(a_grid,d_grid,z_grid,pi_z,n_a,n_d,n_z);

%% Set computational options
vfoptions.do_howard = 1;
vfoptions.n_howard  = 50;
vfoptions.tol       = 1e-6;
vfoptions.verbose   = 0;

muoptions.tol     = 1e-9; 
muoptions.maxiter = 50000;
muoptions.verbose = 0;

% Add these options to structure par
par.vfoptions = vfoptions;
par.muoptions = muoptions;

%% Solve model
% Given params, good initial interval for K/L is [5,6]
% or [0.031363,0.045517] for interest rate r
KL_init = [5.0,6.0]; 
tic
%[p_eq,pol,mu,agg,mom] = solve_model(KL_init,par,grids);
[p_eq,pol,mu,agg,mom] = solve_model_toolkit(KL_init,par,grids);
toc

%% Display results

disp('==================================================================')
disp('PARAMETERS')
fprintf('lambda_hsv : %f \n',par.lambda_hsv)
fprintf('tau_hsv    : %f \n',par.tau_hsv)
disp('------------------------------------')
disp('MOMENTS')
fprintf('Corr(h,z)  : %f \n',mom.corr_h_z)
fprintf('CV(h)      : %f \n',mom.cv.hours)
fprintf('Hours      : %f \n',agg.HH)
fprintf('K/Y        : %f \n',agg.KK/agg.YY)
fprintf('w*L/Y      : %f \n',p_eq.w*agg.LL/agg.YY)
fprintf('I/Y        : %f \n',agg.II/agg.YY)
disp('------------------------------------')
disp('CV')
fprintf('CV(Hours)   : %f \n',mom.cv.hours)
fprintf('CV(Earnings): %f \n',mom.cv.earnings)
fprintf('CV(Income)  : %f \n',mom.cv.income)
fprintf('CV(Wealth)  : %f \n',mom.cv.wealth)
disp('------------------------------------')
disp('GINI')
fprintf('Gini(Hours)   : %f \n',mom.gini.hours)
fprintf('Gini(Earnings): %f \n',mom.gini.earnings)
fprintf('Gini(Income)  : %f \n',mom.gini.income)
fprintf('Gini(Wealth)  : %f \n',mom.gini.wealth)
disp('------------------------------------')
disp('CORR')
fprintf('corr(Hours,z) : %f \n',mom.corr_h_z)
fprintf('corr(Wealth,z): %f \n',mom.corr_k_z)
disp('------------------------------------')
disp('SHARES EARNINGS')
fprintf('q1 earnings: %f \n',mom.shares.earnings(1))
fprintf('q2 earnings: %f \n',mom.shares.earnings(2))
fprintf('q3 earnings: %f \n',mom.shares.earnings(3))
fprintf('q4 earnings: %f \n',mom.shares.earnings(4))
fprintf('q5 earnings: %f \n',mom.shares.earnings(5))
disp('------------------------------------')
disp('SHARES WEALTH')
fprintf('q1 wealth: %f \n',mom.shares.wealth(1))
fprintf('q2 wealth: %f \n',mom.shares.wealth(2))
fprintf('q3 wealth: %f \n',mom.shares.wealth(3))
fprintf('q4 wealth: %f \n',mom.shares.wealth(4))
fprintf('q5 wealth: %f \n',mom.shares.wealth(5))

%% Plots

pol_ap = pol.pol_ap;

figure
plot(a_grid,a_grid,'--','LineWidth',2)
hold on
plot(a_grid,pol_ap(:,1),'LineWidth',2)
hold on
plot(a_grid,pol_ap(:,round(n_z/2)),'LineWidth',2)
hold on
plot(a_grid,pol_ap(:,n_z),'LineWidth',2)
xlabel('assets')
ylabel('aprime')
title('Policy function a''(a,z) ' )

figure
plot(a_grid,sum(mu,2),'LineWidth',2)
xlabel('assets')
ylabel('density')
title('Distribution ' )


