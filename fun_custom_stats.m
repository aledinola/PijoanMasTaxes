function custom_stats = fun_custom_stats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_z,heteroagentoptions,vfoptions,simoptions)
%FUN_CUSTOM_STATS Custom scalar moments for the VFI Toolkit GE routine.
%
% Inputs
%   V,Policy               : value function and policy arrays from the toolkit
%   StationaryDist         : stationary distribution over assets and shocks
%   Params,FnsToEvaluate   : parameter structure and moment functions
%   n_d,n_a,n_z            : grid sizes for labor, assets, and productivity
%   d_grid,a_grid,z_gridvals : labor, asset, and productivity grids
%   pi_z                   : productivity transition matrix
%   heteroagentoptions     : toolkit general-equilibrium options
%   vfoptions,simoptions   : toolkit VFI and simulation options
%
% Output
%   custom_stats           : structure with corr_h_z and cv_hours

custom_stats = struct();

z_mat = repmat(z_gridvals',n_a,1);
pol_d = d_grid(squeeze(Policy(1,:,:)));
custom_stats.corr_h_z = fun_corr(pol_d,z_mat,StationaryDist);

AllStats = EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy, ...
    FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions);

custom_stats.cv_hours = AllStats.H.StdDeviation / AllStats.H.Mean;

end %end function
