%% Pijoan-Mas (2006)
% This program replicates most of the quantitative results of
% Josep Pijoan-Mas, 2006. "Precautionary Savings or Working Longer Hours?"
% Specifically, it reproduces Table 1, Table 2 and Figure 1.
%
% The routine HeteroAgentStationaryEqm_Case1 chooses the parameters beta,
% sigma, lambda, nu to match 4 calibration targets and finds the general
% equilibrium, all in one.

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
vfoptions.tolerance     = 1e-9;
vfoptions.maxiter       = 500;
vfoptions.howards       = 80; 
vfoptions.maxhowards    = 500;
vfoptions.howardsgreedy = 0;
vfoptions.howardssparse = 1;
vfoptions.gridinterplayer = 1;
vfoptions.ngridinterp     = 15;
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
tic;
[~, Policy] = ValueFnIter_Case1( ...
    n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, ...
    ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
time_vfi = toc;

fprintf('Time VFI (seconds): %.2f\n\n', time_vfi);

% --- Policy functions in levels
PolicyValues = PolicyInd2Val_Case1(Policy, n_d, n_a, n_z, d_grid, a_grid, vfoptions);

% --- Stationary distribution
StatDist = StationaryDist_Case1(Policy, n_d, n_a, n_z, pi_z, simoptions);

% --- Additional functions to evaluate
FnsToEvaluate.K        = @(d, aprime, a, z) a;    % assets
FnsToEvaluate.L        = @(d, aprime, a, z) z*d;  % labor in effic units
FnsToEvaluate.H        = @(d, aprime, a, z) d;    % hours worked
FnsToEvaluate.earnings = @(d, aprime, a, z, w) w*z*d; % labor earnings
FnsToEvaluate.income   = @(d, aprime, a, z, r, w) w*z*d + r*a; % income

simoptions.nquantiles = 5;
AllStats = EvalFnOnAgentDist_AllStats_Case1( ...
    StatDist, Policy, FnsToEvaluate, Params, [], ...
    n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);

% Gini coefficients
ginic.wealth   = AllStats.K.Gini;
ginic.hours    = AllStats.H.Gini;
ginic.income   = AllStats.income.Gini;
ginic.earnings = AllStats.earnings.Gini;

% Quintile shares of earnings, income, and wealth
shares.earnings = AllStats.earnings.LorenzCurve([20, 40, 60, 80, 100]) ...
                - [0; AllStats.earnings.LorenzCurve([20, 40, 60, 80])];
shares.income   = AllStats.income.LorenzCurve([20, 40, 60, 80, 100]) ...
                - [0; AllStats.income.LorenzCurve([20, 40, 60, 80])];
shares.wealth   = AllStats.K.LorenzCurve([20, 40, 60, 80, 100]) ...
                - [0; AllStats.K.LorenzCurve([20, 40, 60, 80])];

%% Extra statistics

% Unpack policy functions
StatDist = gather(StatDist);

pol_d  = gather(reshape(PolicyValues(1, :, :), [n_a, n_z])); % d(a,z)
pol_ap = gather(reshape(PolicyValues(2, :, :), [n_a, n_z])); % a'(a,z)

% Average hours by quintile
ave_hours = AllStats.H.QuantileMeans;

%% Correlation statistics: can use either user-written function fun_corr or 
% toolkit function EvalFnOnAgentDist_CrossSectionCorr_InfHorz
% Pijoan-Mas computes corr(h,eps) and corr(a,eps) 
% Correlation between hours and productivity shocks
% Correlation between wealth and productivity shocks

% --- Toolkit command for correlation
FnsToEvaluateCorr.hours = @(d, aprime, a, z) d;
FnsToEvaluateCorr.productivity = @(d, aprime, a, z) z;
FnsToEvaluateCorr.wealth = @(d, aprime, a, z) a;
Corr=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StatDist,Policy,FnsToEvaluateCorr, Params,[], n_d, n_a, n_z, d_grid, a_grid, z_grid,simoptions);

corr_h_z = Corr.hours.CorrelationWith.productivity;
corr_a_z = Corr.wealth.CorrelationWith.productivity;

% --- My own function for correlation
z_mat     = repmat(z_grid', n_a, 1);
corr_h_z2 = fun_corr(pol_d, z_mat, StatDist);
corr_a_z2 = fun_corr(repmat(a_grid,[1,n_z]), z_mat, StatDist);

if abs(corr_h_z-corr_h_z2)>1e-6
    warning('corr_h_z')
end
if abs(corr_a_z-corr_a_z2)>1e-6
    warning('corr_a_z')
end

%% Aggregate moments and dispersion

agg.KK = AllStats.K.Mean;
agg.LL = AllStats.L.Mean;
agg.HH = AllStats.H.Mean;
agg.YY = agg.KK^(1 - Params.theta) * agg.LL^Params.theta;
agg.II = Params.delta * agg.KK;

% Coefficients of variation
cv.hours    = AllStats.H.StdDeviation        / AllStats.H.Mean;
cv.earnings = AllStats.earnings.StdDeviation / AllStats.earnings.Mean;
cv.income   = AllStats.income.StdDeviation   / AllStats.income.Mean;
cv.wealth   = AllStats.K.StdDeviation        / AllStats.K.Mean;

% Policy for consumption
pol_c = Model_cons(pol_d, pol_ap, a_grid, z_grid', Params.K_to_L, Params.theta, Params.delta);

%% Display results on screen

disp('==================================================================');
disp('PARAMETERS');
fprintf('sigma   (Coeff of risk aversion)       : %f \n', Params.sigma);
fprintf('nu      (Curvature labor utility)      : %f \n', Params.nu);
fprintf('lambda  (Weight of labor in disutil)   : %f \n', Params.lambda);
fprintf('beta    (Discount factor)              : %f \n', Params.beta);
fprintf('theta   (Labor share in Cobb-Douglas)  : %f \n', Params.theta);
fprintf('delta   (Capital depreciation rate)    : %f \n', Params.delta);
disp('------------------------------------');
disp('GENERAL EQUILIBRIUM PRICES');
fprintf('K_to_L                                  : %f \n', Params.K_to_L);
fprintf('r                                       : %f \n', Params.r);
fprintf('w                                       : %f \n', Params.w);
disp('------------------------------------');
disp('MOMENTS');
fprintf('Corr(h,z)                               : %f \n', corr_h_z);
fprintf('CV(h)                                   : %f \n', cv.hours);
fprintf('Hours                                   : %f \n', agg.HH);
fprintf('K/Y                                     : %f \n', agg.KK / agg.YY);
fprintf('w*L/Y                                   : %f \n', Params.w * agg.LL / agg.YY);
fprintf('I/Y                                     : %f \n', agg.II / agg.YY);
disp('------------------------------------');
disp('CV');
fprintf('CV(Hours)                               : %f \n', cv.hours);
fprintf('CV(Earnings)                            : %f \n', cv.earnings);
fprintf('CV(Income)                              : %f \n', cv.income);
fprintf('CV(Wealth)                              : %f \n', cv.wealth);
disp('------------------------------------');
disp('GINI');
fprintf('Gini(Hours)                             : %f \n', ginic.hours);
fprintf('Gini(Earnings)                          : %f \n', ginic.earnings);
fprintf('Gini(Income)                            : %f \n', ginic.income);
fprintf('Gini(Wealth)                            : %f \n', ginic.wealth);
disp('------------------------------------');
disp('SHARES EARNINGS');
fprintf('q1 earnings                             : %f \n', shares.earnings(1));
fprintf('q2 earnings                             : %f \n', shares.earnings(2));
fprintf('q3 earnings                             : %f \n', shares.earnings(3));
fprintf('q4 earnings                             : %f \n', shares.earnings(4));
fprintf('q5 earnings                             : %f \n', shares.earnings(5));
disp('------------------------------------');
disp('SHARES WEALTH');
fprintf('q1 wealth                               : %f \n', shares.wealth(1));
fprintf('q2 wealth                               : %f \n', shares.wealth(2));
fprintf('q3 wealth                               : %f \n', shares.wealth(3));
fprintf('q4 wealth                               : %f \n', shares.wealth(4));
fprintf('q5 wealth                               : %f \n', shares.wealth(5));
disp('------------------------------------');
disp('SHARES WEALTH');
fprintf('Corr(h,z) hours and product             : %f \n', corr_h_z);
fprintf('Corr(a,z) wealth and product            : %f \n', corr_a_z);


%% Replicate Table 1 of Pijoan-Mas (2006)  (LaTeX)

if do_save == 1
    fid = fopen(fullfile(results_dir, 'table1.tex'), 'w');

    % LaTeX table preamble
    fprintf(fid, '\\begin{table}[!htbp]\n\\centering\n');
    fprintf(fid, '\\caption{Calibration targets and model parameters}\n');
    fprintf(fid, '\\begin{tabular}{llll}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Parameter & Description & Target & Value \\\\\n');
    fprintf(fid, '\\midrule\n');

    % Rows
    fprintf(fid, '$\\sigma$  & Coeff. risk aversion   & corr(h,\\varepsilon)= %.3f & %.3f \\\\\n', ...
        corr_h_z,      Params.sigma);
    fprintf(fid, '$\\nu$     & Inverse elast. leisure & cv(h)= %.3f              & %.3f \\\\\n', ...
        cv.hours,      Params.nu);
    fprintf(fid, '$\\lambda$ & Weight of leisure      & H= %.3f                  & %.3f \\\\\n', ...
        agg.HH,        Params.lambda);
    fprintf(fid, '$\\beta$   & Discount factor        & K/Y= %.3f                & %.3f \\\\\n', ...
        agg.KK/agg.YY, Params.beta);
    fprintf(fid, '$\\theta$  & Labor share            & wL/Y= %.3f               & %.3f \\\\\n', ...
        Params.w*agg.LL/agg.YY, Params.theta);
    fprintf(fid, '$\\delta$  & Capital depreciation   & I/Y= %.3f                & %.3f \\\\\n', ...
        agg.II/agg.YY, Params.delta);

    % Close table
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');

    fclose(fid);
end

%% Replicate Table 2 of Pijoan-Mas (2006)  (LaTeX)

if do_save == 1
    fid = fopen(fullfile(results_dir, 'table2.tex'), 'w');

    % LaTeX preamble
    fprintf(fid, '\\begin{table}[!htbp]\n\\centering\n');
    fprintf(fid, '\\caption{Distributional statistics}\n');
    fprintf(fid, '\\begin{tabular}{llcccccc}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Variable & $cv$ & Gini & $q_1$ & $q_2$ & $q_3$ & $q_4$ & $q_5$ \\\\\n');
    fprintf(fid, '\\midrule\n');

    % Hours
    fprintf(fid, '\\textbf{Hours} &   &  &   &  &  &  &  \\\\\n');
    fprintf(fid, 'Model $E_0$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', ...
        cv.hours, ginic.hours, ave_hours);

    % Earnings
    fprintf(fid, '\\addlinespace\n');
    fprintf(fid, '\\textbf{Earnings} &   &  &   &  &  &  &  \\\\\n');
    fprintf(fid, 'Model $E_0$ & %.2f & %.2f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\\n', ...
        cv.earnings, ginic.earnings, 100*shares.earnings);

    % Wealth
    fprintf(fid, '\\addlinespace\n');
    fprintf(fid, '\\textbf{Wealth} &   &  &   &  &  &  &  \\\\\n');
    fprintf(fid, 'Model $E_0$ & %.2f & %.2f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\\n', ...
        cv.wealth, ginic.wealth, 100*shares.wealth);

    % Close table
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');

    % Footnote
    fprintf(fid, '\\\\[3ex]\n');
    fprintf(fid, ['\\raggedright\\footnotesize{\\textit{Notes.} $cv$ refers to coefficient of variation. ' ...
                  '$q_1, \\dots, q_5$ refer, for earnings and wealth, to the share held by all people in ' ...
                  'the corresponding quintile with respect to the total. However, for hours it is the ' ...
                  'average number of hours worked by people in the corresponding quintile.}\\\\\n']);
    fprintf(fid, '\\normalsize\n');
    fprintf(fid, '\\end{table}\n');

    fclose(fid);
end

%% Plots

% Distribution of assets
figure;
plot(a_grid, sum(StatDist, 2), 'LineWidth', 2);
xlabel('Assets');
ylabel('Density');
title('Distribution of assets');

if do_save == 1
    print(fullfile(results_dir, 'dist_assets'), '-dpng');
end

% Policy function for assets a'(a,z)
figure;
plot(a_grid, a_grid, '--', 'LineWidth', 2); hold on;
plot(a_grid, pol_ap(:, 1),            'LineWidth', 2);
plot(a_grid, pol_ap(:, round(n_z/2)), 'LineWidth', 2);
plot(a_grid, pol_ap(:, n_z),          'LineWidth', 2);
xlabel('Assets');
ylabel('a''');
title('Policy function a''(a,z)');

if do_save == 1
    print(fullfile(results_dir, 'pol_assets'), '-dpng');
end

% Policy function for hours d(a,z)
figure;
plot(a_grid, pol_d(:, 1),            'LineWidth', 2); hold on;
plot(a_grid, pol_d(:, round(n_z/2)), 'LineWidth', 2);
plot(a_grid, pol_d(:, n_z),          'LineWidth', 2);
legend('z_1', 'z_{mid}', 'z_{high}', 'Location', 'best');
xlabel('Assets');
ylabel('Hours worked');
title('Policy function d(a,z)');

if do_save == 1
    print(fullfile(results_dir, 'pol_hours'), '-dpng');
end

% Policy function for consumption c(a,z)
figure;
plot(a_grid, pol_c(:, 1),            'LineWidth', 2); hold on;
plot(a_grid, pol_c(:, round(n_z/2)), 'LineWidth', 2);
plot(a_grid, pol_c(:, n_z),          'LineWidth', 2);
legend('z_1', 'z_{mid}', 'z_{high}', 'Location', 'best');
xlabel('Assets');
ylabel('Consumption');
title('Policy function c(a,z)');

if do_save == 1
    print(fullfile(results_dir, 'pol_cons'), '-dpng');
end
