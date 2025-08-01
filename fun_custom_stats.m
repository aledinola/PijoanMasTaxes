function custom_stats = fun_custom_stats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_z,heteroagentoptions,vfoptions,simoptions)
% CustomStats=Aiyagari1994_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_z,heteroagentoptions,vfoptions,simoptions)
% % CustomStats: output a structure with custom model stats by fields with the names you want
% % The inputs to CustomModelStats() are fixed and not something the user can
% % change. They differ based on InfHorz of FHorz, and differ if using PType.
% % Note: Most of them are familiar, only thing that may confuse is
% % 'z_gridvals' which is a joint-grid rather than whatever the user set up.
% CustomStats=struct();
% % Note: CustomStats must be scalar-valued

custom_stats = struct();

% corr_h_z
% cv_h

% --- TARGET 1 
% Correlation between hours and productivity shock
z_mat = repmat(z_gridvals',n_a,1);
pol_d  = d_grid(squeeze(Policy(1,:,:))); % d(a,z)
custom_stats.corr_h_z = fun_corr(pol_d,z_mat,StationaryDist);

AllStats=EvalFnOnAgentDist_AllStats_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions);

% --- TARGET 2
% Coefficient of variation of hours of work
custom_stats.cv_hours = AllStats.H.StdDeviation/AllStats.H.Mean;

end %end function