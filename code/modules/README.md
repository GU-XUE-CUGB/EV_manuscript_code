# Analysis modules

| Public code | Module | Entry |
| --- | --- | --- |
| 1 | Code01_EV_reconstruction | `run.py` |
| 2 | Code02_smoothing_cluster_bootstrap | `run.R`; `run_proxy_correlations.R --current-inputs` |
| 3 | Code03_sampling_density_rarefaction | `run.py` |
| 4A | Code04_measurement_uncertainty | `run_analysis.py` |
| 4B | Code05_chronology_uncertainty | `run_analysis.py` |
| 4C | Code06_joint_uncertainty | `run_analysis.py` |
| 5 | Code08_EV_proxy_biota_Pearson | `Supplementary_Code08.R --current-inputs` |
| 5 | Code09_mean_state_biota_Pearson | `Supplementary_Code09.R` |
| 6 | Root `code/` | `code6_Johansen_rate_screening.R`; `code6_ADF_Johansen_VECM.R` |
| 6 plotting | Code10_ADF_Johansen_VECM | `Plot_Fig4.R` |
| 7 | Code07_latitude | `Code07_latitude_analysis.py` |

Use each module's supplied inputs. Component READMEs give commands and dependencies. Plot-only scripts read saved numerical summaries; numerical entry points calculate results.
