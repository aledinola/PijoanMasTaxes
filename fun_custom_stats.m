function custom_stats = fun_custom_stats(~,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,~,~,~,simoptions)
%FUN_CUSTOM_STATS Custom scalar calibration moments for Pijoan-Mas (2006).
%
% Inputs
%   V,Policy               : value function and policy arrays from the toolkit
%   StationaryDist         : stationary distribution over assets and shocks
%   Params,FnsToEvaluate   : parameter structure and moment functions
%   n_d,n_a,n_z            : grid sizes for labor, assets, and productivity
%   d_grid,a_grid,z_gridvals : labor, asset, and productivity grids
%   pi_z                   : productivity transition matrix
%   caliboptions           : toolkit calibration options
%   vfoptions,simoptions   : toolkit VFI and simulation options
%
% Output
%   custom_stats           : structure with six scalar calibration moments

custom_stats = struct();

simoptions_custom = simoptions;
simoptions_custom.whichstats = [1;0;1;0;0;0;0]; % Mean and StdDeviation only

AllStats = EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy, ...
    FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions_custom);

FnsToEvaluateCorr.hours        = @(d, aprime, a, z) d;
FnsToEvaluateCorr.productivity = @(d, aprime, a, z) z;
Corr = EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz( ...
    StationaryDist,Policy,FnsToEvaluateCorr,Params,[], ...
    n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions);

K = AllStats.K.Mean;
L = AllStats.L.Mean;
H = AllStats.H.Mean;
Y = K^(1 - Params.theta) * L^Params.theta;
I = Params.delta * K;
w = fun_w_from_r(Params.r, Params.theta, Params.delta);

custom_stats.corr_h_z = gather(Corr.hours.CorrelationWith.productivity);
custom_stats.cv_h     = gather(AllStats.H.StdDeviation / H);
custom_stats.H        = gather(H);
custom_stats.K_to_Y   = gather(K / Y);
custom_stats.wL_to_Y  = gather(w * L / Y);
custom_stats.I_to_Y   = gather(I / Y);

end %end function
