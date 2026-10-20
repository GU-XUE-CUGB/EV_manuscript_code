# Code 4B — chronological uncertainty

Run `python run_analysis.py` for the configured 10,000-simulation analyses. Required packages: NumPy, pandas and openpyxl. Inputs are in `01_data_revised/`, scripts in `02_code_revised/`, outputs in `03_results/`; `package_config.json` specifies the drivers.

The record-supported scenario includes 24 records / 4,230 measurements: ODP1209, U1403 and Seymour use independent uniform shifts of ±87, ±37 and ±40 kyr. A separate 1,000-simulation ODP1210 comparison shares the ODP1209 Shatsky shift.

Conservative sensitivity includes 38 active records / 8,843 measurements in six independent groups. ODP1209 and ODP1210 share the Shatsky shift. ODP1262 and ODP1267 use separate uniform ±40-kyr shifts. The specified ODP690 post-boundary segment (original age <66.016 Ma) uses a shared uniform −200 to 0 kyr shift; other rows remain fixed. The 8,843-measurement count includes all measurements in active records.

Measurement values remain fixed. Each simulation reassigns bins from perturbed ages, recalculates record-bin IQRs and reapplies n≥2 and IQR>0 eligibility. Common absolute-timescale systematics are excluded. Master SeedSequence entropy is 20260824; group streams follow the order in the source. Saved pilot and reference outputs support comparison checks; the group-only diagnostic is supplied from the pilot.

Plot command: `Rscript 02_code_revised/plot_FigS7_match_user_FigS6_Arial_20261009.R`. It uses saved summaries and display curves in `reference_figures/` and writes to `04_figures/`. Interval boundaries are not smoothed. The plot requires macOS Quartz, Arial and the R packages listed in the script.
