# Code 4A — measurement uncertainty

Run `python run_analysis.py` for the configured 10,000-simulation analyses. Required packages: NumPy, pandas and openpyxl.

`01_data_revised/` contains retained measurements and parameter tables. `02_code_revised/` contains the analysis scripts and helpers. Relative paths and formal drivers are specified in `package_config.json`; outputs go to `03_results/`. Saved reference and pilot files are used for comparison checks.

Documented propagation includes 11 records and 375 measurements. Conservative sensitivity includes 53 records and 10,015 measurements: the documented component, 19 Tier 2 Gaussian records (4,414 measurements) and 23 Tier 2 bounded records (5,226 measurements). Records without an executable observation-level assignment or quantitative model remain fixed.

Reported 2SD and 2σ fields use σ = abs(error)/2; source SD fields retain their SD; specified ± bounds use uniform draws. ODP1049 δ18O uses its species-specific quadratic reconstructions. The specified Tier 2 δ18O conversion uses −5 °C/‰.

Each simulation starts from raw measurements, applies the supplied uncertainty parameters and restandardizes each whole record with sample SD (`ddof=1`). Ages and the observed eligible/nonzero record-bin mask remain fixed. The supplied observation mappings and source tables are retained as inputs.

Plot command: `Rscript 02_code_revised/plot_FigS6_final.R`. It reads saved bin summaries from `reference_figures/` and writes to `04_figures/`. Display curves use `npreg::ss(spar=0.2)` fitted to bin-level medians; interval boundaries are not smoothed. Plot packages are listed in the script.
