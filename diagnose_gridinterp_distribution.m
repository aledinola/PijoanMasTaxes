%% Diagnose grid interpolation policy and stationary distribution
% This script is for debugging only. It compares the standard coarse-grid
% solution with vfoptions.gridinterplayer=1 for the Pijoan-Mas replication.

clear, clc, format long g

addpath(genpath(fullfile('..','VFIToolkit-matlab')))

n_a = 600;
n_d = 51;
n_z = 7;

Params.beta   = 0.945;
Params.sigma  = 1.458;
Params.nu     = 2.833;
Params.lambda = 0.856;
Params.theta  = 0.64;
Params.delta  = 0.083;

Targets.K_to_Y = 3.00;
Params.K_to_L = Targets.K_to_Y^(1 / Params.theta);
[Params.r, Params.w] = fun_prices(Params.K_to_L, Params.theta, Params.delta);

a_curve = 3.0;
a_min   = 0;
a_max   = 50;
a_grid  = a_min + (a_max - a_min) * (linspace(0, 1, n_a)'.^a_curve);
d_grid  = linspace(0.001, 0.999, n_d)';

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

DiscountFactorParamNames = {'beta'};
ReturnFn = @(d, aprime, a, z, r, sigma, lambda, nu, theta, delta) ...
    Model_ReturnFn(d, aprime, a, z, r, sigma, lambda, nu, theta, delta);

basevf = struct();
basevf.lowmemory     = 0;
basevf.verbose       = 0;
basevf.tolerance     = 1e-9;
basevf.maxiter       = 500;
basevf.howards       = 80;
basevf.maxhowards    = 0;
basevf.howardsgreedy = 0;
basevf.howardssparse = 0;
basevf.ngridinterp   = 15;

basesim = struct();
basesim.tolerance = 1e-9;
basesim.ngridinterp = basevf.ngridinterp;

results0 = run_case(0, basevf, basesim, Params, DiscountFactorParamNames, ...
    ReturnFn, n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z);
results1 = run_case(1, basevf, basesim, Params, DiscountFactorParamNames, ...
    ReturnFn, n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z);

fprintf('\n=== Cross-case policy comparison ===\n');
fprintf('max |hours GI - noGI|: %.12g\n', max(abs(results1.pol_d(:) - results0.pol_d(:))));
fprintf('max |aprime GI - noGI|: %.12g\n', max(abs(results1.pol_ap(:) - results0.pol_ap(:))));
fprintf('mean |hours GI - noGI|: %.12g\n', mean(abs(results1.pol_d(:) - results0.pol_d(:))));
fprintf('mean |aprime GI - noGI|: %.12g\n', mean(abs(results1.pol_ap(:) - results0.pol_ap(:))));

function out = run_case(use_gi, basevf, basesim, Params, DiscountFactorParamNames, ...
    ReturnFn, n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z)

vfoptions = basevf;
vfoptions.gridinterplayer = use_gi;
simoptions = basesim;
simoptions.gridinterplayer = use_gi;

fprintf('\n=== gridinterplayer = %d ===\n', use_gi);
[~, Policy] = ValueFnIter_InfHorz(n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, ...
    ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
PolicyValues = PolicyInd2Val_InfHorz(Policy, n_d, n_a, n_z, d_grid, a_grid, vfoptions);
StatDist = StationaryDist_InfHorz(Policy, n_d, n_a, n_z, pi_z, simoptions);

FnsToEvaluate.K = @(d, aprime, a, z) a;
FnsToEvaluate.L = @(d, aprime, a, z) z*d;
FnsToEvaluate.H = @(d, aprime, a, z) d;
FnsToEvaluate.earnings = @(d, aprime, a, z, w) w*z*d;
FnsToEvaluate.income = @(d, aprime, a, z, r, w) w*z*d + r*a;
FnsToEvaluate.C = @(d, aprime, a, z, r, theta, delta) ...
    Model_cons(d, aprime, a, z, r, theta, delta);
simoptions.nquantiles = 5;
AllStats = EvalFnOnAgentDist_AllStats_InfHorz(StatDist, Policy, FnsToEvaluate, ...
    Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);

FnsToEvaluateCorr.hours = @(d, aprime, a, z) d;
FnsToEvaluateCorr.productivity = @(d, aprime, a, z) z;
FnsToEvaluateCorr.wealth = @(d, aprime, a, z) a;
Corr = EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StatDist, Policy, ...
    FnsToEvaluateCorr, Params, [], n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions);

Policy_cpu = gather(Policy);
PolicyValues_cpu = gather(PolicyValues);
StatDist_cpu = gather(StatDist);

out.pol_d = reshape(PolicyValues_cpu(1,:,:), [n_a,n_z]);
out.pol_ap = reshape(PolicyValues_cpu(2,:,:), [n_a,n_z]);
out.StatDist = StatDist_cpu;

dist_mass = sum(StatDist_cpu(:));
mean_a = sum(StatDist_cpu(:) .* repmat(a_grid, [n_z,1]));
mean_ap = sum(StatDist_cpu(:) .* out.pol_ap(:));
fprintf('Policy size: %s\n', mat2str(size(Policy_cpu)));
fprintf('dist mass: %.16f, min: %.12g, max: %.12g, negatives: %d\n', ...
    dist_mass, min(StatDist_cpu(:)), max(StatDist_cpu(:)), sum(StatDist_cpu(:) < -1e-14));
fprintf('mean a from dist: %.12g\n', mean_a);
fprintf('mean aprime under dist: %.12g\n', mean_ap);
Y = AllStats.K.Mean^(1 - Params.theta) * AllStats.L.Mean^Params.theta;
fprintf('moments K/Y: %.12g, wL/Y: %.12g, I/Y: %.12g, CV(H): %.12g, Gini(K): %.12g\n', ...
    AllStats.K.Mean / Y, Params.w * AllStats.L.Mean / Y, ...
    Params.delta * AllStats.K.Mean / Y, ...
    AllStats.H.StdDeviation / AllStats.H.Mean, AllStats.K.Gini);
fprintf('corr(h,z): %.12g, corr(a,z): %.12g\n', ...
    Corr.hours.CorrelationWith.productivity, Corr.wealth.CorrelationWith.productivity);
fprintf('mass at first five a-grid points: %s\n', mat2str(sum(StatDist_cpu(1:5,:),2)', 12));
fprintf('mass at last five a-grid points:  %s\n', mat2str(sum(StatDist_cpu(end-4:end,:),2)', 12));
fprintf('hours range: [%.12g, %.12g]\n', min(out.pol_d(:)), max(out.pol_d(:)));
fprintf('aprime range: [%.12g, %.12g]\n', min(out.pol_ap(:)), max(out.pol_ap(:)));

if use_gi == 1
    raw_d = reshape(Policy_cpu(1,:,:), [n_a,n_z]);
    raw_l1 = reshape(Policy_cpu(2,:,:), [n_a,n_z]);
    raw_l2 = reshape(Policy_cpu(3,:,:), [n_a,n_z]);
    raw_flag = reshape(Policy_cpu(4,:,:), [n_a,n_z]);

    fine_grid = interp1(gpuArray(1:n_a)', gpuArray(a_grid), ...
        linspace(1, n_a, n_a + (n_a - 1) * vfoptions.ngridinterp)');
    fine_grid = gather(fine_grid);
    fine_index = (vfoptions.ngridinterp + 1) * (raw_l1 - 1) + raw_l2;
    decoded_ap = fine_grid(fine_index);
    fprintf('raw d index range: [%g, %g]\n', min(raw_d(:)), max(raw_d(:)));
    fprintf('raw L1 range: [%g, %g], raw L2 range: [%g, %g]\n', ...
        min(raw_l1(:)), max(raw_l1(:)), min(raw_l2(:)), max(raw_l2(:)));
    fprintf('L2 flag counts [1 2 3]: [%d %d %d]\n', ...
        sum(raw_flag(:)==1), sum(raw_flag(:)==2), sum(raw_flag(:)==3));
    fprintf('max |PolicyInd2Val aprime - decoded fine aprime|: %.12g\n', ...
        max(abs(out.pol_ap(:) - decoded_ap(:))));

    dist_usual = stationary_from_policy_gi(raw_l1, raw_l2, raw_flag, pi_z, ...
        vfoptions.ngridinterp, false);
    dist_flag = stationary_from_policy_gi(raw_l1, raw_l2, raw_flag, pi_z, ...
        vfoptions.ngridinterp, true);
    mean_a_usual = sum(dist_usual(:) .* repmat(a_grid, [n_z,1]));
    mean_a_flag = sum(dist_flag(:) .* repmat(a_grid, [n_z,1]));
    fprintf('manual dist, ignoring L2 flag, mean a: %.12g\n', mean_a_usual);
    fprintf('manual dist, using L2 flag, mean a:    %.12g\n', mean_a_flag);
    fprintf('max |toolkit dist - manual usual|: %.12g\n', max(abs(StatDist_cpu(:) - dist_usual(:))));
    fprintf('max |toolkit dist - manual flag|:  %.12g\n', max(abs(StatDist_cpu(:) - dist_flag(:))));
    fprintf('manual flag mass at first five a-grid points: %s\n', ...
        mat2str(sum(dist_flag(1:5,:),2)', 12));
end

end %end function

function dist = stationary_from_policy_gi(raw_l1, raw_l2, raw_flag, pi_z, ngridinterp, use_flag)

[n_a, n_z] = size(raw_l1);
N = n_a * n_z;

prob_upper = (raw_l2 - 1) / (ngridinterp + 1);
if use_flag
    prob_upper(raw_flag == 1) = 0;
    prob_upper(raw_flag == 3) = 1;
end
prob_lower = 1 - prob_upper;

row = zeros(2 * n_z * N, 1);
col = zeros(2 * n_z * N, 1);
val = zeros(2 * n_z * N, 1);
ctr = 0;

for iz = 1:n_z
    current = (1:n_a)' + n_a * (iz - 1);
    lower = raw_l1(:,iz);
    upper = min(lower + 1, n_a);
    for izp = 1:n_z
        ctr = ctr + n_a;
        idx = ctr-n_a+1:ctr;
        row(idx) = lower + n_a * (izp - 1);
        col(idx) = current;
        val(idx) = prob_lower(:,iz) * pi_z(iz,izp);

        ctr = ctr + n_a;
        idx = ctr-n_a+1:ctr;
        row(idx) = upper + n_a * (izp - 1);
        col(idx) = current;
        val(idx) = prob_upper(:,iz) * pi_z(iz,izp);
    end
end

T = sparse(row, col, val, N, N);
dist = zeros(N, 1);
dist(ceil(n_a/2):n_a:N) = 1 / n_z;
dist = dist / sum(dist);

for it = 1:200000
    newdist = T * dist;
    if max(abs(newdist - dist)) < 1e-12
        dist = newdist;
        break
    end
    dist = newdist;
end

dist = reshape(dist / sum(dist), [n_a,n_z]);

end %end function
