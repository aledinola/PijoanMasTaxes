function [Params,pol,StationaryDist,agg,mom] = solve_model_toolkit(par,grids,vfoptions,simoptions,heteroagentoptions)

% p_eq: KL,r,w
% agg: KK,LL,HH,YY
% mom: cv,corr_h_z,corr_k_z,gini,shares

%% Unpack grids

a_grid = grids.a_grid;
d_grid = grids.d_grid;
z_grid = grids.z_grid;
pi_z   = grids.pi_z;
n_a  = grids.n_a;
n_d  = grids.n_d;
n_z  = grids.n_z;

Params.crra   = par.crra;
Params.nu     = par.nu;
Params.lambda = par.lambda;
Params.beta   = par.beta;
Params.theta  = par.theta;
Params.delta  = par.delta;
Params.lambda_hsv = par.lambda_hsv;
Params.tau_hsv    = par.tau_hsv;

% Initial guess for GE parameter
Params.K_to_L = par.K_to_L;

%% Setup toolkit inputs
DiscountFactorParamNames={'beta'};

ReturnFn = @(d,aprime,a,z,K_to_L,crra,lambda,nu,theta,delta,lambda_hsv,tau_hsv) Model_ReturnFn(d,aprime,a,z,K_to_L,crra,lambda,nu,theta,delta,lambda_hsv,tau_hsv);

%% Solve general equil 

% Create functions to be evaluated
FnsToEvaluate.K = @(d,aprime,a,z) a;   % A, assets or capital
FnsToEvaluate.L = @(d,aprime,a,z) z*d; % L, labor in efficiency units
FnsToEvaluate.H = @(d,aprime,a,z) d;   % H, Hours of work

% Now define the functions for the General Equilibrium conditions
    % Should be written as LHS of general eqm eqn minus RHS, so that the closer the value given by the function is to 
    % zero, the closer the general eqm condition is to holding.
GeneralEqmEqns.CapitalMarket = @(K_to_L,K,L) K_to_L-K/L; %The requirement that the interest rate corresponds to the agg capital level
% Inputs can be any parameter, price, or aggregate of the FnsToEvaluate

GEPriceParamNames={'K_to_L'};
% Set initial value for K/L


% Solve for the stationary general equilbirium

fprintf('Calculating price vector corresponding to the stationary general eqm \n')
[p_eqm,~,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

disp(GeneralEqmCondn)

%% Recompute at equil prices
if heteroagentoptions.maxiter>0
    Params.K_to_L = p_eqm.K_to_L; 
else
    Params.K_to_L = par.K_to_L;
end
[Params.r, Params.w] = fun_prices(Params.K_to_L,Params);

fprintf('Calculating various equilibrium objects \n')
tic
[~,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
time_vfi = toc;
disp(' ')
fprintf('Time VFI (seconds): %.2f',time_vfi)
disp(' ')
% V is value function
% Policy is policy function (but as an index of k_grid, not the actual values)

% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_s,d_grid,a_grid); % This will give you the policy in terms of values rather than index

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

% Functions to be evaluated
FnsToEvaluate.K = @(d,aprime,a,z) a; 
FnsToEvaluate.L = @(d,aprime,a,z) z*d;
FnsToEvaluate.H = @(d,aprime,a,z) d;
FnsToEvaluate.hours    = @(d,aprime,a,z,w) d;
FnsToEvaluate.earnings = @(d,aprime,a,z,w) w*z*d;
FnsToEvaluate.income   = @(d,aprime,a,z,r,w) w*z*d+r*a;
FnsToEvaluate.wealth   = @(d,aprime,a,z,r,w) a;

%AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Parameters, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions, EntryExitParamNames, PolicyWhenExiting)
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate,Params, [],n_d, n_a, n_z, d_grid, a_grid,z_grid,simoptions);

% AggVars contains the aggregate values of the 'FnsToEvaluate' (in this model aggregates are equal to the mean expectation over the agent distribution)
% Currently the only FnsToEvaluate is assets, so we get aggregate capital stock

%                                                StationaryDist, PolicyIndexes, FnsToEvaluate, Parameters, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions
AllStats = EvalFnOnAgentDist_AllStats_Case1(StationaryDist, Policy, FnsToEvaluate, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid,simoptions);

%TopWealthShares=100*(1-LorenzCurves.Wealth([80,95,99],1)); % Need the 20,5, and 1 top shares for Tables of Huggett (1996)

% Calculate the wealth gini
gini_wealth   = AllStats.wealth.Gini;
gini_hours    = AllStats.hours.Gini;
gini_income   = AllStats.income.Gini;
gini_earnings = AllStats.earnings.Gini;

% Quintile Shares of Earnings, Income, and Wealth
shares.earnings = AllStats.earnings.LorenzCurve([20,40,60,80,100])-[0;AllStats.earnings.LorenzCurve([20,40,60,80])] ;
shares.income   = AllStats.income.LorenzCurve([20,40,60,80,100])-[0;AllStats.income.LorenzCurve([20,40,60,80])] ;
shares.wealth   = AllStats.wealth.LorenzCurve([20,40,60,80,100])-[0;AllStats.wealth.LorenzCurve([20,40,60,80])] ;

% % Calculate cross-sectional correlations
% FnsToEvaluateExtra.wealth = @(d,aprime,a,z) a;
% FnsToEvaluateExtra.hours  = @(d,aprime,a,z) d;
% FnsToEvaluateExtra.z      = @(d,aprime,a,z) z;
% 
% CrossSectionCorr = EvalFnOnAgentDist_CrossSectionCorr_Case1(StationaryDist, Policy, FnsToEvaluateExtra, Params,[], n_d, n_a, n_z, d_grid, a_grid, z_grid);

%% Prepare Outputs

pol_ind_d  = squeeze(Policy(1,:,:)); % d(a,z)
pol_ind_ap = squeeze(Policy(2,:,:)); % a'(a,z)

pol_d  = d_grid(pol_ind_d);
pol_ap = a_grid(pol_ind_ap);

pol = pack_into_struct(pol_d,pol_ap);

% Aggregate moments
theta  = par.theta;
delta  = par.delta;
agg.KK = AggVars.K.Mean;
agg.LL = AggVars.L.Mean;
agg.HH = AggVars.H.Mean;
agg.YY = agg.KK^(1-theta)*agg.LL^theta;
agg.II = delta*agg.KK;

% Coefficient of variation for hours, earnings, income and wealth
mom.cv.hours    = AllStats.H.StdDeviation/AllStats.H.Mean;
mom.cv.earnings = AllStats.earnings.StdDeviation/AllStats.earnings.Mean;
mom.cv.income   = AllStats.income.StdDeviation/AllStats.income.Mean;
mom.cv.wealth   = AllStats.wealth.StdDeviation/AllStats.wealth.Mean;

% Gini coeff for hours, earnings, income and wealth
mom.gini.hours     = gini_hours;
mom.gini.earnings  = gini_earnings;
mom.gini.income    = gini_income;
mom.gini.wealth    = gini_wealth;

% Cross-sectional correlations
% Order: (k,h,z)
mom.corr_h_z  =  NaN;%CrossSectionCorr.hours.z;
mom.corr_k_z  =  NaN;%  CrossSectionCorr.wealth.z;

mom.shares = shares;

end %end function solve_model