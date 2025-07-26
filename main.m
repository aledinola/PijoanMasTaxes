%% Pijoan-Mas (2006)
clear
clc
close all
format long g
% This is the folder where the VFI toolkit files are saved
mypath = 'C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab';
addpath(genpath(mypath))

%% Set computational options
do_pijoan = 1;   % If 1, load shocks from Pijoan-Mas files, otherwise discretize
n_a       = 300; % No. grid points for assets
n_d       = 51;  % No. grid points for labor supply

% --- Value functions options
vfoptions=struct(); 
vfoptions.lowmemory     = 1;
vfoptions.verbose       = 0;
vfoptions.tolerance     = 1e-6;
vfoptions.maxiter       = 500;
vfoptions.howards       = 80; 
vfoptions.maxhowards    = 500;
vfoptions.howardsgreedy = 0;
vfoptions.gridinterplayer = 1;
vfoptions.ngridinterp     = 15;

% Distribution options
simoptions=struct(); % Use default options for solving for stationary distribution
simoptions.gridinterplayer = vfoptions.gridinterplayer;
simoptions.ngridinterp     = vfoptions.ngridinterp;

% Heteroagentoptions
heteroagentoptions = struct();
heteroagentoptions.verbose=1; % verbose means that you want it to give you feedback on what is going on
heteroagentoptions.toleranceGEprices=1e-4; % default is 1e-4
heteroagentoptions.toleranceGEcondns=1e-4; % default is 1e-4
heteroagentoptions.fminalgo = 1;
heteroagentoptions.maxiter = 0;

%% Set economic parameters

par.crra   = 1.458; % Coeff of risk aversion
par.nu     = 2.833; % Curvature labor utility
par.lambda = 0.856; % Weight of labor in util
par.beta   = 0.945; % Discount factor
par.theta  = 0.64;  % Labor share in Cobb-Douglas
par.delta  = 0.083; % Capital depreciation rate

% Initial guess for GE parameter/price
par.K_to_L = 5.5345;

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

a_grid = make_grid(a_min,a_max,n_a,a_curve,1);

%grid for labor
d_grid = linspace(0.001,0.999,n_d)';

if do_pijoan==1
    disp('Pijoan-Mas discretization of z with 7 points')
    fprintf('n_a = %d, n_d = %d \n',n_a,n_d)
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
    disp('AR(1) discretization of z with n_z points')
    n_z = 7;
    [pi_z,z_grid_log] = markovapprox(rho_z,sig_z,0.0,3.0,n_z,0);
    z_grid     = exp(z_grid_log);

end

% Pack grids into a struct
grids = pack_into_struct(a_grid,d_grid,z_grid,pi_z,n_a,n_d,n_z);


%% Solve model
%tic
[Params,pol,mu,agg,mom] = solve_model_toolkit(par,grids,vfoptions,simoptions,heteroagentoptions);
%toc

%% Display results

disp('==================================================================')
disp('PARAMETERS')
fprintf('crra (Coeff of risk aversion)       : %f \n',par.crra) 
fprintf('nu (Curvature labor utility)        : %f \n',par.nu)
fprintf('lambda (Weight of labor in disutil) : %f \n',par.lambda)
fprintf('beta (Discount factor)              : %f \n',par.beta) 
fprintf('theta (Labor share in Cobb-Douglas) : %f \n',par.theta) 
fprintf('delta (Capital depreciation rate)   : %f \n',par.delta) 
fprintf('lambda_hsv : %f \n',par.lambda_hsv)
fprintf('tau_hsv    : %f \n',par.tau_hsv)
disp('------------------------------------')
disp('GENERAL EQUILIBRIUM PRICES')
fprintf('K_to_L : %f \n',Params.K_to_L)
fprintf('r      : %f \n',Params.r)
fprintf('w      : %f \n',Params.w)
disp('------------------------------------')
disp('MOMENTS')
fprintf('Corr(h,z)  : %f \n',mom.corr_h_z)
fprintf('CV(h)      : %f \n',mom.cv.hours)
fprintf('Hours      : %f \n',agg.HH)
fprintf('K/Y        : %f \n',agg.KK/agg.YY)
fprintf('w*L/Y      : %f \n',Params.w*agg.LL/agg.YY)
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
pol_d  = pol.pol_d;

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

figure
plot(a_grid,pol_d(:,1),'LineWidth',2)
hold on
plot(a_grid,pol_d(:,round(n_z/2)),'LineWidth',2)
hold on
plot(a_grid,pol_d(:,n_z),'LineWidth',2)
legend('z_1','z_4','z_7')
xlabel('Assets')
ylabel('Hours worked')
title('Policy function d(a,z) ' )


