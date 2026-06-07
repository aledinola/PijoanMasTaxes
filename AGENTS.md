# Repository Instructions

This repo replicates Pijoan-Mas (2006), "Precautionary Savings or Working
Longer Hours?" The paper is the source of truth for model equations and
calibration targets. The replication document
[results/PijoanMas2006_main.tex](results/PijoanMas2006_main.tex) summarizes the
paper equations and compares the MATLAB outputs with the original tables and
figure, so check it before changing model equations, targets, or reported
objects.

The code uses functions from the VFI Toolkit, available on this computer at
`C:\Users\aledi\OneDrive\Documents\GitHub\VFIToolkit-matlab` and online at
[https://github.com/vfitoolkit](https://github.com/vfitoolkit). It is also a
compact test case for VFI Toolkit workflows:

1. Infinite-horizon general-equilibrium solving, with and without the
   GridInterpolationLayer. The `gridinterplayer` switch is a field of both
   `vfoptions` and `simoptions`; the previous issue with `gridinterplayer=1`
   was fixed in the toolkit, so treat this as a supported option rather than a
   known problem.
2. General-equilibrium calibration. The `do_calib` flag in [main.m](main.m)
   controls whether the script solves once at fixed parameters or calibrates
   parameters in general equilibrium by minimizing the distance between model
   moments and data moments while also solving GE conditions.

## MATLAB Structure

- [main.m](main.m) is the entry point. Start there for paths, options,
  parameters, grids, shocks, VFI Toolkit inputs, GE conditions, solve vs
  calibration flags, statistics, tables, and plots.
- [Model_ReturnFn.m](Model_ReturnFn.m) defines household period utility.
- [Model_cons.m](Model_cons.m) defines consumption from the household budget
  constraint.
- [fun_prices.m](fun_prices.m) maps capital-labor ratios to prices.
- [fun_w_from_r.m](fun_w_from_r.m) gives the wage, and optionally the
  capital-labor ratio, implied by an interest rate.
- [fun_calibration_moments.m](fun_calibration_moments.m) defines the six scalar
  calibration moments shared by calibration and final reporting.
- [fun_custom_stats.m](fun_custom_stats.m) is the `CalibrateBIHAModel` callback
  that computes toolkit statistics and calls `fun_calibration_moments`.
- [fun_corr.m](fun_corr.m) is a simple weighted-correlation helper.

## Main Workflow

In [main.m](main.m), use `do_GE` to choose general-equilibrium solving and
`do_calib` to choose calibration. When `do_calib == 0` and `do_GE == 1`, the
script solves the stationary equilibrium at the current parameter values. When
`do_calib == 1`, it calls `CalibrateBIHAModel` to update selected parameters
and the equilibrium interest rate jointly.

For clarity, all calibration targets are kept together as
`TargetMoments.CustomModelStats`: `corr_h_z`, `cv_h`, `H`, `K_to_Y`, `wL_to_Y`,
and `I_to_Y`. `fun_custom_stats.m` uses a local `simoptions.whichstats` setting
for calibration speed; the final reporting call to
`EvalFnOnAgentDist_AllStats_InfHorz` in `main.m` should keep full statistics
available for tables and plots. When `do_save == 1`, `main.m` writes
`results/calibration_summary.txt` with parameter values, data targets, final
model moments, weights, and the calibration objective value when available.

The `gridinterplayer` option should be set consistently through
`vfoptions.gridinterplayer` and `simoptions.gridinterplayer`. Use
`gridinterplayer=1` for the finer-grid interpolation workflow and
`gridinterplayer=0` for the coarse-grid baseline.

## Replication Document

[results/PijoanMas2006_main.tex](results/PijoanMas2006_main.tex) is organized
to reproduce Table 1, Table 2, and Figure 1 from Pijoan-Mas (2006). For each
object, the document shows the original object from Pijoan-Mas and the
corresponding MATLAB replication output.

Do not re-extract or reread
[pijoan_mas_paper/PijoanMas_RED_2006.pdf](pijoan_mas_paper/PijoanMas_RED_2006.pdf)
unless existing extracted artifacts are missing or insufficient. To save time
and tokens, store any PDF-reading artifacts, such as extracted text snippets,
page-range text files, or rendered pages, inside `pijoan_mas_paper`.

## Toolkit Boundary

The VFI Toolkit is stored outside this repo at
`C:\Users\aledi\OneDrive\Documents\GitHub\VFIToolkit-matlab`; the upstream
project is [https://github.com/vfitoolkit](https://github.com/vfitoolkit).
Inspect the local toolkit checkout when behavior depends on toolkit internals.
The most relevant toolkit folders for this replication are `ValueFnIter`,
`StationaryDist`, `HeterogeneousAgent`, `EvaluateFnOnAgentDist`, `PolicyInd2Val`,
`SubCodes`, `SimulateTimeSeries`, and `Estimation`.

For infinite-horizon GE solving, start from
`HeterogeneousAgent/InfHorz/HeteroAgentStationaryEqm_InfHorz.m`. For GE
calibration, start from `Estimation/Calibration/CalibrateBIHAModel.m`; its
joint and nested objective functions are
`Estimation/ObjectiveFn/CalibrateBIHAModel_Joint_objectivefn.m` and
`Estimation/ObjectiveFn/CalibrateBIHAModel_Nested_objectivefn.m`. The
GridInterpolationLayer implementation is documented at
[robertdkirkby/GridInterpolationLayer](https://github.com/robertdkirkby/GridInterpolationLayer)
and is also available locally at
`C:\Users\aledi\OneDrive\Documents\GitHub\GridInterpolationLayer`.

Codex may freely modify files in this replication repo. Before modifying any
toolkit file, Codex must ask permission and create or switch to a toolkit-repo
branch named `ale`, so changes stay separate from Robert Kirkby's official repo
history.
