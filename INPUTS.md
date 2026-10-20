# Input files

| Module | Numerical inputs |
| --- | --- |
| Code 1 | Compiled measurement inventory in `Code01_EV_reconstruction/inputs/`; reconstruction comparisons in `reference_results/` |
| Code 2 | Record-IQR trajectories and bin means in `inputs/`; fits/envelopes in `reference_results/` and `plot_snapshot_inputs/`; current EV–proxy plotting tables in `reference_results/current_dataset/` |
| Code 3 | Retained measurements, record-bin IQRs and site mapping in `inputs/`; summaries in `reference_results/` and `plot_snapshot_inputs/` |
| Code 4A–C | Raw measurements in `01_data_revised/`; parameter/source workbooks and pilot acceptance references in `03_results/EV_uncertainty_inventory/`; reference arrays in `reference_outputs/`; plot summaries in `reference_figures/` |
| Code 5 | Current environmental–biological pairs and environmental tables in `Code08_EV_proxy_biota_Pearson/inputs/current_*.csv`; mean-state inputs in `Code09_mean_state_biota_Pearson/inputs/` |
| Code 6 | Extinction paired datasets embedded in `code6_ADF_Johansen_VECM.R`; three-rate screening uses the current Code 5 pair table; Fig. 4 statistics in `code/results/code6_current_dataset/` |
| Code 7 | Record-bin IQRs, latitude/site mappings and the S1 workbook in `Code07_latitude/inputs/`; numerical plot inputs in `output/`; biological/event artwork in `assets/` |

Module paths are relative to `code/modules/`, except for explicit `code/` paths.

Record-bin IQRs require at least 2 measurements; positive IQRs contribute to EV means. Measurement-only simulations retain the observed eligibility mask; chronology and joint simulations reassign bins and reapply eligibility.

Keep the supplied source tables and parameter workbooks with the scripts. Uncertainty folders named `pilot`, `source_inputs` or `Chronology_shared_age_model_grouping_audit` contain required inputs or acceptance references.

The Code 5 root entry selects `--current-inputs`. Direct variability-component commands must also use that flag: `Rscript code/modules/Code08_EV_proxy_biota_Pearson/Supplementary_Code08.R --current-inputs` and `Rscript code/modules/Code08_EV_proxy_biota_Pearson/Supplementary_pCO2_variability.R --current-inputs`.
