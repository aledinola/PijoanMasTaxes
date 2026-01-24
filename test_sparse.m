%% Test option howardssparse on model with d variable and gridinterplayer=1

clear,clc,close all,format long g

%% Paths and saving options

% Folder where the VFI toolkit files are saved
mypath = 'C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab';
addpath(genpath(mypath))

% Flag for saving output (tables and figures)
do_save     = 0;              % Set to 1 to save LaTeX tables and figures
results_dir = 'results';      % Folder for saved output

if do_save == 1 && ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Set computational options

do_GE     = 0;
do_pijoan = 1;   % If 1, load shocks from Pijoan-Mas files, otherwise discretize
n_a       = 600; % No. grid points for assets
n_d       = 51;  % No. grid points for labor supply

% --- Value functions options
vfoptions=struct(); 
vfoptions.lowmemory     = 0;
vfoptions.verbose       = 0;
vfoptions.tolerance     = 1e-10;
vfoptions.maxiter       = 500;
vfoptions.howards       = 80; 
vfoptions.maxhowards    = 500;
vfoptions.howardsgreedy = 0;
vfoptions.howardssparse = 0;
vfoptions.gridinterplayer = 1;
vfoptions.ngridinterp     = 15;
vfoptions.maxaprimediff   = 30;
%vfoptions.divideandconquer = 0;

% --- Distribution / simulation options
simoptions                  = struct(); % Use default options for stationary distribution
simoptions.tolerance        = 1e-9;
simoptions.gridinterplayer  = vfoptions.gridinterplayer;
simoptions.ngridinterp      = vfoptions.ngridinterp;
%simoptions.inheritanceasset = 0;

% --- Heterogeneous-agent GE options
heteroagentoptions                          = struct();
heteroagentoptions.verbose                  = 1;      % 1 = print progress
heteroagentoptions.toleranceGEprices        = 1e-4;   % default 1e-4
heteroagentoptions.toleranceGEcondns        = 1e-4;   % default 1e-4
heteroagentoptions.fminalgo                 = 1;      % 0=fzero, 1=fminsearch, 8=lsqnonlin
heteroagentoptions.maxiter                  = 1000;

%% Set economic parameters

% Initial guess for GE parameter/price
Params.K_to_L = 5.5649;

% --- Preference parameters
Params.beta   = 0.945; % Discount factor
Params.sigma  = 1.458; % Coefficient of relative risk aversion
Params.nu     = 2.833; % Curvature of labor disutility
Params.lambda = 0.856; % Weight of labor in utility

% --- Technology parameters
Params.theta  = 0.64;  % Labor share in Cobb–Douglas
Params.delta  = 0.083; % Capital depreciation rate

% --- Parameters for AR(1) labor productivity z (not used if do_pijoan = 1)
Params.rho_z = 0.92;
Params.sig_z = 0.21;

%% Set grids and shocks

% Grid for assets
a_curve = 3.0;
a_min   = 0;
a_max   = 50;
a_grid  = a_min + (a_max - a_min) * (linspace(0, 1, n_a)'.^a_curve);

% Grid for labor
d_grid  = linspace(0.001, 0.999, n_d)';

% Productivity shocks
if do_pijoan == 1
    disp('Pijoan-Mas discretization of z with 7 points');
    fprintf('n_a = %d, n_d = %d \n', n_a, n_d);

    z_grid_log = [-1.385493 -0.923662 -0.461831 0.000000 0.461831 0.923662 1.385493]';
    z_grid     = exp(z_grid_log);

    pi_z = [ ...
        0.746464  0.252884  0.000652  0.000000  0.000000  0.000000  0.000000
        0.046088  0.761085  0.192512  0.000314  0.000000  0.000000  0.000000
        0.000028  0.069422  0.788612  0.141793  0.000145  0.000000  0.000000
        0.000000  0.000065  0.100953  0.797965  0.100953  0.000065  0.000000
        0.000000  0.000000  0.000145  0.141793  0.788612  0.069422  0.000028
        0.000000  0.000000  0.000000  0.000314  0.192512  0.761085  0.046088
        0.000000  0.000000  0.000000  0.000000  0.000652  0.252884  0.746464];

    pi_z = pi_z ./ sum(pi_z, 2);
    n_z  = length(z_grid);

else
    disp('AR(1) discretization of z with n_z points');
    n_z = 7;
    [pi_z, z_grid_log] = markovapprox(Params.rho_z, Params.sig_z, 0.0, 3.0, n_z, 0);
    z_grid              = exp(z_grid_log);
end

%% Setup toolkit inputs

DiscountFactorParamNames = {'beta'};

% --- Model payoff function
ReturnFn = @(d, aprime, a, z, K_to_L, sigma, lambda, nu, theta, delta) ...
    Model_ReturnFn(d, aprime, a, z, K_to_L, sigma, lambda, nu, theta, delta);

% --- Functions to evaluate on the distribution
FnsToEvaluate.K = @(d, aprime, a, z) a;    % Assets / capital
FnsToEvaluate.L = @(d, aprime, a, z) z .* d; % Labor in efficiency units
FnsToEvaluate.H = @(d, aprime, a, z) d;    % Hours of work

% --- General equilibrium conditions
% Conditions are written as LHS - RHS = 0 at equilibrium.
GeneralEqmEqns.CapitalMarket = @(K_to_L, K, L) K_to_L - K ./ L;

% beta calibrated to match K/Y = 3.0
GeneralEqmEqns.target_beta   = @(K, L, theta) K ./ (K.^(1 - theta) .* L.^theta) - 3.0;

% lambda calibrated to match H = 0.33
GeneralEqmEqns.target_lambda = @(H) H - 0.33;

% sigma (CRRA) calibrated to match corr(hours, z) = 0.02
GeneralEqmEqns.target_sigma  = @(corr_h_z) corr_h_z - 0.02;

% nu calibrated to match coeff. of variation of hours = 0.22
GeneralEqmEqns.target_nu     = @(cv_hours) cv_hours - 0.22;

GEPriceParamNames                        = {'K_to_L', 'beta', 'lambda', 'sigma', 'nu'};
heteroagentoptions.constrain0to1         = {'beta'};
heteroagentoptions.multiGEweights        = [1, 1, 1, 1, 1];

heteroagentoptions.CustomModelStats = @(V, Policy, StationaryDist, Params, ...
    FnsToEvaluate, n_d, n_a, n_z, d_grid, a_grid, z_gridvals, pi_z, ...
    heteroagentoptions, vfoptions, simoptions) ...
    fun_custom_stats(V, Policy, StationaryDist, Params, FnsToEvaluate, ...
                     n_d, n_a, n_z, d_grid, a_grid, z_gridvals, pi_z, ...
                     heteroagentoptions, vfoptions, simoptions);

%% Solve for the stationary general equilibrium (optional)

if do_GE == 1
    fprintf('Calculating the stationary general equilibrium...\n');

    [p_eqm, ~, GeneralEqmCondn] = HeteroAgentStationaryEqm_Case1( ...
        n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ...
        ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, ...
        DiscountFactorParamNames, [], [], [], ...
        GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);

    disp('GeneralEqmCondn:');
    disp(GeneralEqmCondn);

    % Update GE parameters in structure Params
    for ii = 1:numel(GEPriceParamNames)
        name_ii       = GEPriceParamNames{ii};
        Params.(name_ii) = p_eqm.(name_ii);
    end
end

%% Recompute at equilibrium prices

[Params.r, Params.w] = fun_prices(Params.K_to_L, Params.theta, Params.delta);

fprintf('Calculating value functions, policies, and stationary distribution...\n');

% --- Value function iteration
vfoptions.howardssparse=0;
tic;
[V, Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
time_vfi = toc;
fprintf('Time VFI (seconds) with howardssparse=0: %.2f\n\n', time_vfi);

vfoptions.howardssparse=1;
tic;
[V1, Policy1]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
time_vfi = toc;
fprintf('Time VFI (seconds) with howardssparse=1: %.2f\n\n', time_vfi);

err_P = max(abs(Policy(:)-Policy1(:)))
err_V = max(abs(V(:)-V1(:)))

% --- Policy functions in levels
PolicyValues = PolicyInd2Val_Case1(Policy, n_d, n_a, n_z, d_grid, a_grid, vfoptions);
pol_d  = gather(reshape(PolicyValues(1, :, :), [n_a, n_z])); % d(a,z)
pol_ap = gather(reshape(PolicyValues(2, :, :), [n_a, n_z])); % a'(a,z)

%% Plots
% Policy function for assets a'(a,z)
figure;
plot(a_grid, a_grid, '--', 'LineWidth', 2); hold on;
plot(a_grid, pol_ap(:, 1),            'LineWidth', 2);
plot(a_grid, pol_ap(:, round(n_z/2)), 'LineWidth', 2);
plot(a_grid, pol_ap(:, n_z),          'LineWidth', 2);
xlabel('Assets');
ylabel('a''');
title('Policy function a''(a,z)');

% Policy function for hours d(a,z)
figure;
plot(a_grid, pol_d(:, 1),            'LineWidth', 2); hold on;
plot(a_grid, pol_d(:, round(n_z/2)), 'LineWidth', 2);
plot(a_grid, pol_d(:, n_z),          'LineWidth', 2);
legend('z_1', 'z_{mid}', 'z_{high}', 'Location', 'best');
xlabel('Assets');
ylabel('Hours worked');
title('Policy function d(a,z)');