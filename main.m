%% Pijoan-Mas (2006)
% This program replicates most of the quantitative results of
% Josep Pijoan-Mas, 2006. "Precautionary Savings or Working Longer Hours?"
% Specifically, it reproduces Table 1, Table 2 and Figure 1.
%
% The script can either solve the model at fixed parameter values or 
% calibrate selected parameters in general equilibrium. The GE search is 
% written in terms of r, which is equivalent to searching over K/L through 
% the firm's FOC.

clear,clc,close all,format long g

%% Paths and saving options

% Folder where the VFI toolkit files are saved
mypath = fullfile('..','VFIToolkit-matlab');
addpath(genpath(mypath))

% Flag for saving output (tables, figures, and calibration summary)
do_save     = 1;              % Set to 1 to save outputs
results_dir = 'results';      % Folder for saved output

if do_save == 1 && ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Set computational options

do_calib  = 0;   % 0=solve once, 1=calibrate six parameters and r in GE
do_GE     = 0;   % 0=solve at fixed r, 1=solve GE over r
do_pijoan = 1;   % If 1, load shocks from Pijoan-Mas files
n_a       = 600; % No. grid points for assets
n_d       = 51;  % No. grid points for labor supply

% --- Value functions options
vfoptions=struct(); 
vfoptions.lowmemory     = 0;    % If 1, less memory intensive
vfoptions.verbose       = 0;
vfoptions.tolerance     = 1e-9; % VFI tolerance
vfoptions.maxiter       = 500;  % VFI max number of iterations
vfoptions.howards       = 80;   % Number of iterations for Howard 
vfoptions.maxhowards    = 0;
vfoptions.howardsgreedy = 0;
vfoptions.howardssparse = 0;
vfoptions.gridinterplayer = 0;  % 0=a' on coarse grid, 1=a' on finer grid
vfoptions.ngridinterp     = 15; % number of extra points between coarse-grid points
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
heteroagentoptions.maxiter                  = 100;

%% Set economic parameters

% --- Preference parameters: initial guesses for calibration, or fixed values
%     when do_calib == 0. Trailing comments report the paper values.
Params.beta   = 0.96; % paper: 0.945, discount factor
Params.sigma  = 1.3;  % paper: 1.458, coefficient of relative risk aversion
Params.nu     = 2.5;  % paper: 2.833, curvature of labor disutility
Params.lambda = 0.5;  % paper: 0.856, weight of labor in utility

% --- Technology parameters
Params.theta  = 0.64;  % Labor share in Cobb-Douglas
Params.delta  = 0.083; % Capital depreciation rate

% Paper targets used for calibration and data comparisons.
Targets.corr_h_z = 0.02;
Targets.cv_h     = 0.22;
Targets.H        = 0.33;
Targets.K_to_Y   = 3.00;
Targets.wL_to_Y  = 0.64;
Targets.I_to_Y   = 0.25;
CalibMomentNames = {'corr_h_z','cv_h','H','K_to_Y','wL_to_Y','I_to_Y'};
CalibParamNames = {'beta','sigma','nu','lambda','theta','delta'};
CalibWeights = ones(numel(CalibMomentNames),1);

% Fixed-price runs start from the paper target K/Y=3. In GE and calibration
% runs, K/Y=(K/L)^theta is used only to construct the initial r.
Params.K_to_L = Targets.K_to_Y^(1 / Params.theta);
[Params.r, Params.w] = fun_prices(Params.K_to_L, Params.theta, Params.delta);

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

% Productivity shocks: use the paper discretization by default, otherwise
% approximate the AR(1) process with Tauchen's method.
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
    Tauchen_q = 3.0;
    %[pi_z, z_grid_log] = markovapprox(Params.rho_z, Params.sig_z, 0.0, 3.0, n_z, 0);
    [z_grid_log,pi_z] = discretizeAR1_Tauchen(0.0,Params.rho_z,Params.sig_z,n_z,Tauchen_q);
    z_grid              = exp(z_grid_log);
end

%% Setup toolkit inputs

DiscountFactorParamNames = {'beta'};

% --- Model payoff function and GE price object
% The GE variable is r. This is equivalent to searching over K/L because the
% firm's FOC maps r one-to-one into K/L.
ReturnFn = @(d, aprime, a, z, r, sigma, lambda, nu, theta, delta) ...
    Model_ReturnFn(d, aprime, a, z, r, sigma, lambda, nu, theta, delta);

% --- Minimal functions needed for GE conditions and calibration targets
FnsToEvaluate.K = @(d, aprime, a, z) a;   % Assets / capital
FnsToEvaluate.L = @(d, aprime, a, z) z*d; % Labor in efficiency units
FnsToEvaluate.H = @(d, aprime, a, z) d;   % Hours of work

% --- General equilibrium condition
% The condition is guessed r minus firm-implied r. In stationarity, aggregate
% current assets and next-period assets coincide.
GeneralEqmEqns.CapitalMarket = @(r,K,L,theta,delta) r-((1-theta)*(K / L)^(-theta)-delta);

GEPriceParamNames = {'r'};

%% Calibrate or solve for the stationary general equilibrium

if do_calib == 1
    fprintf('Calibrating the stationary general equilibrium...\n');

    ParametrizeParamsFn = [];

    TargetMoments = struct();
    TargetMoments.CustomModelStats.corr_h_z = Targets.corr_h_z;
    TargetMoments.CustomModelStats.cv_h     = Targets.cv_h;
    TargetMoments.CustomModelStats.H        = Targets.H;
    TargetMoments.CustomModelStats.K_to_Y   = Targets.K_to_Y;
    TargetMoments.CustomModelStats.wL_to_Y  = Targets.wL_to_Y;
    TargetMoments.CustomModelStats.I_to_Y   = Targets.I_to_Y;

    caliboptions = struct();
    caliboptions.verbose = 1;
    caliboptions.CustomModelStats = @fun_custom_stats;
    caliboptions.jointoptimization = 1;
    caliboptions.fminalgo = 8;
    caliboptions.metric = 'sum_squared';
    caliboptions.weights = CalibWeights;
    caliboptions.constrainAtoB = CalibParamNames(:);
    caliboptions.constrainAtoBlimits.beta   = [0.800, 0.995];
    caliboptions.constrainAtoBlimits.sigma  = [1.010, 5.000];
    caliboptions.constrainAtoBlimits.nu     = [1.010, 10.00];
    caliboptions.constrainAtoBlimits.lambda = [0.010, 5.000];
    caliboptions.constrainAtoBlimits.theta  = [0.100, 0.950];
    caliboptions.constrainAtoBlimits.delta  = [0.001, 0.200];

    calib_heteroagentoptions = heteroagentoptions;
    calib_heteroagentoptions.constrainAtoB = GEPriceParamNames(:);
    calib_heteroagentoptions.constrainAtoBlimits.r = [0.001, 0.150];

    tStart = tic;
    [CalibParams, calibsummary] = CalibrateBIHAModel( ...
        CalibParamNames, TargetMoments, ...
        n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, ...
        ReturnFn, Params, DiscountFactorParamNames, ParametrizeParamsFn, ...
        GEPriceParamNames, FnsToEvaluate, GeneralEqmEqns, ...
        calib_heteroagentoptions, caliboptions, vfoptions, simoptions);
    time_calib = toc(tStart);
    fprintf('CalibrateBIHAModel: %.4f seconds\n', time_calib);

    calib_param_names = fieldnames(CalibParams);
    for ii = 1:numel(calib_param_names)
        name_ii = calib_param_names{ii};
        Params.(name_ii) = CalibParams.(name_ii);
    end

    fprintf('Calibration objective value: %g\n', calibsummary.objvalue);

end

if do_calib == 0 && do_GE == 1
    fprintf('Calculating the stationary general equilibrium...\n');

    [p_eqm, ~, GeneralEqmCondn] = HeteroAgentStationaryEqm_InfHorz( ...
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

[Params.w, Params.K_to_L] = fun_w_from_r(Params.r, Params.theta, Params.delta);

fprintf('Calculating value functions, policies, and stationary distribution...\n');

% --- Value function iteration
tic;
[~, Policy] = ValueFnIter_InfHorz(n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, ...
    ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
time_vfi = toc;

fprintf('Time VFI (seconds): %.2f\n\n', time_vfi);

% --- Policy functions in levels
PolicyValues = PolicyInd2Val_InfHorz(Policy, n_d, n_a, n_z, d_grid, a_grid, vfoptions);

% --- Stationary distribution
StatDist = StationaryDist_InfHorz(Policy, n_d, n_a, n_z, pi_z, simoptions);

% --- Full set of functions for final reporting
FnsToEvaluate.K        = @(d, aprime, a, z) a;    % assets
FnsToEvaluate.L        = @(d, aprime, a, z) z*d;  % labor in efficiency units
FnsToEvaluate.H        = @(d, aprime, a, z) d;    % hours worked
FnsToEvaluate.earnings = @(d, aprime, a, z, w) w*z*d; % labor earnings
FnsToEvaluate.income   = @(d, aprime, a, z, r, w) w*z*d + r*a; % income
FnsToEvaluate.C        = @(d, aprime, a, z, r, theta, delta) Model_cons(d, aprime, a, z, r, theta, delta); % consumption

simoptions.nquantiles = 5;
AllStats = EvalFnOnAgentDist_AllStats_InfHorz(StatDist, Policy, FnsToEvaluate,...
    Params,[],n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);

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

% Unpack policy functions for hours and a'
StatDist = gather(StatDist);
pol_d  = gather(reshape(PolicyValues(1, :, :), [n_a, n_z])); % d(a,z)
pol_ap = gather(reshape(PolicyValues(2, :, :), [n_a, n_z])); % a'(a,z)

% Policy for consumption
ValOnGrid=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy,FnsToEvaluate,Params,...
    [], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);
pol_c = ValOnGrid.C; % c(a,z)

if ~isequal(size(pol_c),[n_a,n_z])
    error('Policy for consumption computed by ValuesOnGrid is not correct')
end

% Average hours by quintile
ave_hours = AllStats.H.QuantileMeans;

%% Correlation statistics

% --- Toolkit cross-sectional correlations for reporting
FnsToEvaluateCorr.hours        = @(d, aprime, a, z) d;
FnsToEvaluateCorr.productivity = @(d, aprime, a, z) z;
FnsToEvaluateCorr.wealth       = @(d, aprime, a, z) a;
Corr=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StatDist,Policy,FnsToEvaluateCorr, Params,[], n_d, n_a, n_z, d_grid, a_grid, z_grid,simoptions);

ModelMoments = fun_calibration_moments(AllStats,Corr,Params);

corr_h_z = ModelMoments.corr_h_z;
corr_a_z = gather(Corr.wealth.CorrelationWith.productivity);

%% Aggregate moments and dispersion

agg.KK = AllStats.K.Mean;
agg.LL = AllStats.L.Mean;
agg.HH = AllStats.H.Mean;
agg.YY = agg.KK^(1 - Params.theta) * agg.LL^Params.theta;
agg.II = Params.delta * agg.KK;
GE_resid = GeneralEqmEqns.CapitalMarket(Params.r, agg.KK, agg.LL, ...
    Params.theta, Params.delta);

% Coefficients of variation
cv.hours    = ModelMoments.cv_h;
cv.earnings = AllStats.earnings.StdDeviation / AllStats.earnings.Mean;
cv.income   = AllStats.income.StdDeviation   / AllStats.income.Mean;
cv.wealth   = AllStats.K.StdDeviation        / AllStats.K.Mean;

%% Final calibration target moments

if do_save == 1
    fid = fopen(fullfile(results_dir, 'calibration_summary.txt'), 'w');

    fprintf(fid, 'Parameter values\n');
    fprintf(fid, '----------------\n');
    summary_param_names = [CalibParamNames, GEPriceParamNames];
    for ii = 1:numel(summary_param_names)
        name_ii = summary_param_names{ii};
        fprintf(fid, '%-12s %.12g\n', name_ii, Params.(name_ii));
    end

    fprintf(fid, '\nCalibration targets\n');
    fprintf(fid, '-------------------\n');
    fprintf(fid, '%-12s %16s %16s %16s\n', 'moment', 'data', 'model', 'weight');
    for ii = 1:numel(CalibMomentNames)
        name_ii = CalibMomentNames{ii};
        fprintf(fid, '%-12s %16.12g %16.12g %16.12g\n', ...
            name_ii, Targets.(name_ii), ModelMoments.(name_ii), CalibWeights(ii));
    end

    fprintf(fid, '\nCalibration objective value\n');
    fprintf(fid, '---------------------------\n');
    if do_calib == 1
        fprintf(fid, '%.12g\n', calibsummary.objvalue);
    else
        fprintf(fid, 'not computed (do_calib=0)\n');
    end

    fclose(fid);
end

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
fprintf('Hours                                   : %f \n', ModelMoments.H);
fprintf('K/Y                                     : %f \n', ModelMoments.K_to_Y);
fprintf('w*L/Y                                   : %f \n', ModelMoments.wL_to_Y);
fprintf('I/Y                                     : %f \n', ModelMoments.I_to_Y);
fprintf('GE capital-market residual              : %e \n', GE_resid);
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
disp('CORRELATIONS');
fprintf('Corr(h,z) hours and product             : %f \n', corr_h_z);
fprintf('Corr(a,z) wealth and product            : %f \n', corr_a_z);


%% Replicate Table 1 of Pijoan-Mas (2006)  (LaTeX)

if do_save == 1
    fid = fopen(fullfile(results_dir, 'table1.tex'), 'w');

    % LaTeX table preamble
    fprintf(fid, '\\begin{table}[!htbp]\n\\centering\n');
    fprintf(fid, '\\caption{Calibration targets and model parameters}\n');
    fprintf(fid, '\\begin{tabular}{lcc}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Parameter & Target & Value \\\\\n');
    fprintf(fid, '\\midrule\n');

    fprintf(fid, '$\\sigma$  & $\\operatorname{corr}(h,\\varepsilon)=%.2f$ & %.3f \\\\\n', ...
        ModelMoments.corr_h_z, Params.sigma);
    fprintf(fid, '$\\nu$     & $cv(h)=%.2f$ & %.3f \\\\\n', ...
        ModelMoments.cv_h, Params.nu);
    fprintf(fid, '$\\lambda$ & $H=%.2f$ & %.3f \\\\\n', ...
        ModelMoments.H, Params.lambda);
    fprintf(fid, '$\\beta$   & $K/Y=%.2f$ & %.3f \\\\\n', ...
        ModelMoments.K_to_Y, Params.beta);
    fprintf(fid, '$\\theta$  & $wL/Y=%.2f$ & %.3f \\\\\n', ...
        ModelMoments.wL_to_Y, Params.theta);
    fprintf(fid, '$\\delta$  & $I/Y=%.2f$ & %.3f \\\\\n', ...
        ModelMoments.I_to_Y, Params.delta);

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
    fprintf(fid, '\\begin{tabular}{lccccccc}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Variable & $cv$ & Gini & $q_1$ & $q_2$ & $q_3$ & $q_4$ & $q_5$ \\\\\n');
    fprintf(fid, '\\midrule\n');

    % Hours
    fprintf(fid, '\\multicolumn{8}{l}{Hours} \\\\\n');
    fprintf(fid, 'Model $E_0$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', ...
        cv.hours, ginic.hours, ave_hours);
    fprintf(fid, 'Data (CPS) & 0.22 & 0.11 & 0.24 & 0.31 & 0.33 & 0.35 & 0.42 \\\\\n');

    % Earnings
    fprintf(fid, '\\addlinespace\n');
    fprintf(fid, '\\multicolumn{8}{l}{Earnings} \\\\\n');
    fprintf(fid, 'Model $E_0$ & %.2f & %.2f & %.1f\\%% & %.1f\\%% & %.1f\\%% & %.1f\\%% & %.1f\\%% \\\\\n', ...
        cv.earnings, ginic.earnings, 100*shares.earnings);
    fprintf(fid, 'Data (CPS) & 0.56 & 0.29 & 7.9\\%% & 13.7\\%% & 18.0\\%% & 23.3\\%% & 37.1\\%% \\\\\n');
    fprintf(fid, 'Data (SCF) & 2.65 & 0.61 & -0.2\\%% & 4.0\\%% & 13.0\\%% & 22.9\\%% & 60.2\\%% \\\\\n');

    % Wealth
    fprintf(fid, '\\addlinespace\n');
    fprintf(fid, '\\multicolumn{8}{l}{Wealth} \\\\\n');
    fprintf(fid, 'Model $E_0$ & %.2f & %.2f & %.1f\\%% & %.1f\\%% & %.1f\\%% & %.1f\\%% & %.1f\\%% \\\\\n', ...
        cv.wealth, ginic.wealth, 100*shares.wealth);
    fprintf(fid, 'Data (SCF) & 6.53 & 0.80 & -0.3\\%% & 1.3\\%% & 5.0\\%% & 12.2\\%% & 81.7\\%% \\\\\n');

    % Close table
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');

    % Footnote
    fprintf(fid, '\\\\[3ex]\n');
    fprintf(fid, ['\\raggedright\\footnotesize{\\textit{Notes.} $cv$ refers to coefficient of variation. ' ...
                  '$q_1, \\dots, q_5$ refer, for earnings and wealth, to the share held by all people in ' ...
                  'the corresponding quintile with respect to the total. However, for hours it is the ' ...
                  'average number of hours worked by people in the corresponding quintile. Statistics ' ...
                  'from SCF correspond to the 1998 wave and are quoted from Budria et al. (2002). ' ...
                  'Statistics from CPS correspond to the 2002 wave.}\\\\\n']);
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
plot(a_grid, pol_d(:, 1), 'k-', 'LineWidth', 1.2); hold on;
plot(a_grid, pol_d(:, 4), 'k--', 'LineWidth', 1.2);
plot(a_grid, pol_d(:, 7), 'k:', 'LineWidth', 1.8);
legend('\epsilon_1', '\epsilon_4', '\epsilon_7', 'Location', 'northeast');
xlabel('assets');
ylabel('hours');
title('Hours');
xlim([0, 50]);
ylim([0, 0.60]);
box on;

if do_save == 1
    print(fullfile(results_dir, 'pol_hours'), '-dpng');
end

% Policy function for consumption c(a,z)
figure;
plot(a_grid, pol_c(:, 1),  'k-',          'LineWidth', 2); hold on;
plot(a_grid, pol_c(:, round(n_z/2)),'k--', 'LineWidth', 2);
plot(a_grid, pol_c(:, n_z), 'k:',      'LineWidth', 2);
legend('z_1', 'z_{mid}', 'z_{high}', 'Location', 'best');
xlabel('Assets');
ylabel('Consumption');
title('Policy function c(a,z)');

if do_save == 1
    print(fullfile(results_dir, 'pol_cons'), '-dpng');
end
