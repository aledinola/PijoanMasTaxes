# GE Calibration Flag For `main.m`

**Summary**
- Add `do_calib` near the existing run flags in `main.m`.
- Default `do_calib = 0` keeps current behavior: solve once in PE or GE depending on `do_GE`.
- When `do_calib = 1`, call `CalibrateBIHAModel` for GE calibration of `beta`, `sigma`, `nu`, `lambda`, `theta`, and `delta`.
- Use six targets: `corr_h_z`, `cv_h`, `H`, `K_to_Y`, `wL_to_Y`, and `I_to_Y`.

**Implementation Changes**
- In `main.m`, activate `Targets.corr_h_z = 0.02` and keep the other five targets already present.
- Translate flat `Targets` into toolkit-facing:
  `TargetMoments.CustomModelStats.corr_h_z`, `cv_h`, `H`, `K_to_Y`, `wL_to_Y`, `I_to_Y`.
- Add:
  `CalibParamNames = {'beta','sigma','nu','lambda','theta','delta'};`
  `ParametrizeParamsFn = [];`
  `caliboptions.CustomModelStats = @fun_custom_stats;`
  `caliboptions.jointoptimization = 1;`
  `caliboptions.fminalgo = 8;`
- Use bounded transforms:
  `beta [0.80, 0.995]`, `sigma [1.01, 5]`, `nu [1.01, 10]`, `lambda [0.01, 5]`, `theta [0.10, 0.95]`, `delta [0.001, 0.20]`.
- Bound joint GE price `r` with `heteroagentoptions.constrainAtoB = {'r'}` and limits `[0.001, 0.15]`.
- After calibration, copy `CalibParams` fields back into `Params`, including joint-optimized `r`, then continue through the existing recompute/reporting block.

**Custom Moments**
- Update `fun_custom_stats.m` to return exactly:
  `corr_h_z`, `cv_h`, `H`, `K_to_Y`, `wL_to_Y`, `I_to_Y`.
- Compute correlations only with toolkit correlation machinery:
  `EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz`.
- Remove `fun_corr` from the calibration path and from the post-solve correlation check in `main.m`; keep toolkit correlation as the single source for `corr_h_z` and `corr_a_z`.
- Compute `cv_h` and `H` from toolkit `AllStats`.
- Compute `K_to_Y`, `wL_to_Y`, and `I_to_Y` from aggregate `K`, `L`, calibrated `theta`, `delta`, and the wage implied by calibrated `r`.

**Reporting**
- Keep Table 1’s “Target” column populated from realized model moments, matching the current MATLAB style.
- Print `calibsummary.objvalue` after calibration.
- After final recompute, print the final GE capital-market residual.

**Test Plan**
- Run default `do_calib = 0` to confirm current single-solve behavior still works.
- Run a reduced-grid calibration smoke test with `do_calib = 1`, `do_save = 0`, and smaller grids.
- Run full calibration only if runtime is reasonable; otherwise report that only the reduced-grid calibration path was verified.

**Assumptions**
- Calibration and estimation are synonyms here.
- `do_calib = 1` always means GE calibration, regardless of `do_GE`.
- Use joint `lsqnonlin`; avoid nested calibration unless toolkit code is separately fixed with permission.
