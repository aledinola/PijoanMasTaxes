# Repository Instructions

This repo replicates Pijoan-Mas (2006), "Precautionary Savings or Working
Longer Hours?" The source of truth is the paper. Its equations are summarized
in [results/PijoanMas2006_main.tex](results/PijoanMas2006_main.tex), so check
that file before changing model equations or calibration targets.

Use this replication to test two VFI Toolkit workflows:

1. Infinite-horizon general-equilibrium solving, with and without linear
   interpolation. The `gridinterplayer` switch is a field of both `vfoptions`
   and `simoptions`. The relevant toolkit entry point is
   `HeterogeneousAgent/InfHorz/HeteroAgentStationaryEqm_InfHorz.m`.
2. General-equilibrium calibration. In this repo, calibration means choosing
   parameter values to minimize the weighted distance between model moments and
   data moments while also solving GE conditions. The relevant toolkit entry
   point is `Estimation/Calibration/CalibrateBIHAModel.m`, but this function
   is not well tested and may contain bugs.

## MATLAB Structure

- [main.m](main.m) is the entry point. Start there for paths, options,
  parameters, grids, shocks, VFI Toolkit inputs, GE conditions, equilibrium
  solving, statistics, tables, and plots.
- [Model_ReturnFn.m](Model_ReturnFn.m) defines household period utility.
- [Model_cons.m](Model_cons.m) defines consumption from the household budget
  constraint.
- [fun_prices.m](fun_prices.m) maps capital-labor ratios to prices.
- [fun_w_from_r.m](fun_w_from_r.m) gives the wage, and optionally the
  capital-labor ratio, implied by an interest rate.
- [fun_custom_stats.m](fun_custom_stats.m) and [fun_corr.m](fun_corr.m) compute
  custom moments used for reporting and calibration-style comparisons.

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

## Branch Workflow

The `main` branch is intentionally left showing the current
`gridinterplayer=1` problem, so VFI Toolkit users can reproduce and inspect the
issue. Use the `ale` branch as the temporary investigation playground for
experiments, diagnostics, and fixes.

## Toolkit Boundary

The VFI Toolkit is stored outside this repo at
`C:\Users\aledi\OneDrive\Documents\GitHub\VFIToolkit-matlab`. Inspect it when
behavior depends on toolkit internals. The most relevant toolkit folders for
this replication are `ValueFnIter`, `StationaryDist`, `HeterogeneousAgent`,
`EvaluateFnOnAgentDist`, `PolicyInd2Val`, `SubCodes`, `SimulateTimeSeries`, and
`Estimation`.

For infinite-horizon GE solving, start from
`HeterogeneousAgent/InfHorz/HeteroAgentStationaryEqm_InfHorz.m`. For GE
calibration, start from `Estimation/Calibration/CalibrateBIHAModel.m`; its
joint and nested objective functions are
`Estimation/ObjectiveFn/CalibrateBIHAModel_Joint_objectivefn.m` and
`Estimation/ObjectiveFn/CalibrateBIHAModel_Nested_objectivefn.m`. Treat this
calibration path as less mature than the standard solver path, and verify it
carefully when using it.

Codex may freely modify files in this replication repo. Before modifying any
toolkit file, Codex must ask permission and create or switch to a toolkit-repo
branch named `ale`, so changes stay separate from Robert Kirkby's official repo
history.

Here there is some detailed info on how grisinterplayer is done in the toolkit: [robertdkirkby/GridInterpolationLayer](https://github.com/robertdkirkby/GridInterpolationLayer),
also available on my machine here: "C:\Users\aledi\OneDrive\Documents\GitHub\GridInterpolationLayer"
