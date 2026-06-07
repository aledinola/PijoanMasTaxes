function custom_stats = fun_custom_stats(~,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,~,~,~,simoptions)
%FUN_CUSTOM_STATS Custom scalar calibration moments for Pijoan-Mas (2006).
%
% Inputs
%   Policy                 : toolkit policy array at the current parameters
%   StationaryDist         : stationary distribution over assets and shocks
%   Params,FnsToEvaluate   : parameter structure and moment functions
%   n_d,n_a,n_z            : grid sizes for labor, assets, and productivity
%   d_grid,a_grid,z_gridvals : labor, asset, and productivity grids
%   simoptions             : toolkit simulation options; copied locally before
%                            restricting the requested AllStats outputs
%
% Output
%   custom_stats           : structure with six scalar calibration moments

simoptions_custom = simoptions;
simoptions_custom.whichstats = [1;0;1;0;0;0;0]; % Mean and StdDeviation only

AllStats = EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy, ...
    FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions_custom);

FnsToEvaluateCorr.hours        = @(d, aprime, a, z) d;
FnsToEvaluateCorr.productivity = @(d, aprime, a, z) z;
Corr = EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz( ...
    StationaryDist,Policy,FnsToEvaluateCorr,Params,[], ...
    n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions);

custom_stats = fun_calibration_moments(AllStats,Corr,Params);

end %end function
