# Auxiliary sources

- `SST/SST_smooth spline_portable.R` reads `SST/inputs/SST.csv` and smooths an already calibrated SST series. Screening, calibration and 45° correction are upstream input steps.
- `nannofossils/variant2/N_smooth spline_predict_portable.R` uses the included Lowery2020 extinction-percentage input. It smooths extinction percentages and predicts on a 0.1-Myr grid; it does not estimate speciation or net diversification.

These scripts write to their local `generated/` directories. The supplied workflow-diagram script is described in `workflow_figures/README.md`.
