# Code 4C — joint measurement and chronological uncertainty

Run `python run_analysis.py` for the configured 10,000-simulation analysis. Required packages: NumPy, pandas and openpyxl. Inputs are in `01_data_revised/`, scripts in `02_code_revised/`, outputs in `03_results/`; `package_config.json` specifies the drivers.

The analysis combines the conservative measurement parameters of Code 4A and the six-group conservative chronology parameters of Code 4B. Common absolute-timescale systematics are excluded.

Master SeedSequence entropy 20260824 spawns independent measurement [0] and chronology [1] streams. The chronology stream spawns [1,0] through [1,5] in the specified group order.

Each simulation begins with raw retained values and ages, applies both perturbations, recomputes whole-record z-scores using sample SD (`ddof=1`), reassigns bins and reapplies n≥2 and IQR>0 eligibility. Saved pilot and single-component summaries are comparison inputs, not substitutes for generating joint trajectories.

Plot command: `Rscript 02_code_revised/plot_FigS8_joint_conservative_measurement_chronology_uncertainty.R`. It reads the supplied formal summary from `reference_figures/` and writes to `04_figures/`. Display curves use `npreg::ss(spar=0.2)` fitted to 49 bin medians; interval boundaries remain unsmoothed. Plot packages are listed in the script.
