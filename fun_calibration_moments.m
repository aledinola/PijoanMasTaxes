function ModelMoments = fun_calibration_moments(AllStats,Corr,Params)
%FUN_CALIBRATION_MOMENTS Six scalar moments used for GE calibration.
%
% Inputs
%   AllStats               : toolkit statistics structure for K, L, and H
%   Corr                   : toolkit cross-sectional correlation structure
%   Params                 : parameter structure containing r, theta, and delta
%
% Output
%   ModelMoments           : structure with calibration moment scalars

K = AllStats.K.Mean;
L = AllStats.L.Mean;
H = AllStats.H.Mean;
Y = K^(1 - Params.theta) * L^Params.theta;

% Return CPU scalars because calibration and reporting assemble plain vectors.
ModelMoments.corr_h_z = gather(Corr.hours.CorrelationWith.productivity);
ModelMoments.cv_h     = gather(AllStats.H.StdDeviation / H);
ModelMoments.H        = gather(H);
ModelMoments.K_to_Y   = gather(K / Y);
ModelMoments.wL_to_Y  = gather(fun_w_from_r(Params.r, Params.theta, Params.delta) * L / Y);
ModelMoments.I_to_Y   = gather(Params.delta * K / Y);

end %end function
