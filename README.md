# K–Pg environmental variability: Supplementary Code

Analysis code for *Escalating environmental instability exacerbated the end-Cretaceous marine ecosystem collapse*.

## Contents

| Code | Analysis | Related figures / datasets |
| --- | --- | --- |
| 1 | Standardization, binning, IQR and EV aggregation | Datasets S1–S6 |
| 2 | Smoothing, background, record-cluster bootstrap and EV–proxy associations | Fig. 2; S3, S10–S11 |
| 3 | Sampling density and rarefaction | Fig. S4–S5 |
| 4A | Measurement uncertainty | Fig. S6 |
| 4B | Chronological uncertainty | Fig. S7 |
| 4C | Joint measurement and chronological uncertainty | Fig. S8 |
| 5 | Environmental variability / mean state and biological Pearson correlations | Fig. 3; S12–S13; Dataset S9 |
| 6 | ADF, Johansen, VECM and directional tests | Fig. 4; Dataset S10 |
| 7 | Latitudinal EV, site-balanced rarefaction and leave-one-site-out analysis | Fig. 5; S14–S15 |

Numbered scripts are in `code/`; functions and inputs are in `code/modules/`. Keep the directory structure. Internal module folder numbers are mapped to the public Code 1–7 order in `code/modules/README.md`.

## Requirements

Python and R are required. Packages are listed in `requirements.txt` and `requirements-R.txt`; saved `R_sessionInfo.txt` files record R environments. Arial/Quartz plotting scripts require Arial and macOS graphics support, as specified in their component READMEs.

Install Python dependencies with `python3 -m pip install -r requirements.txt`. Install the listed R packages separately.

## Analysis commands

Run from the repository root. These commands execute analyses and may be computationally intensive. Each component uses its supplied inputs; do not automatically replace another component's input snapshot.

```sh
# Code 1: compiled-inventory reconstruction
python3 code/modules/Code01_EV_reconstruction/run.py

# Code 2: smoothing, 5,000 bootstrap resamples and EV–proxy associations
Rscript "code/code2_smoothing spline.R" --cluster-bootstrap --proxy-correlations

# Code 3: 10,000 rarefaction iterations
python3 "code/code3_bootstrap resampling_all data.py"

# Code 4A–C: 10,000 simulations for each configured scenario
python3 "code/code4_grouped resampling_all data.py"
python3 code/code4B_chronological_uncertainty.py
python3 code/code4C_joint_measurement_chronology.py

# Code 5: biological Pearson analyses and mean environmental states
Rscript "code/code5_Pearsen analysis.R"

# Code 6: three-rate screening and the two extinction VECMs
Rscript code/code6_Johansen_rate_screening.R

# Code 7: site-balanced rarefaction and figure generation
Rscript "code/code7_bootstrap resampling_mid-high latitude.R"
```

Code 6 screening sources `code6_ADF_Johansen_VECM.R`, which first executes the two extinction models; a separate extinction run is not needed. Code 4B runs its record-supported scenario before the conservative scenario, which uses the first scenario's generated draw archive. Code 7 uses `python3`, or the interpreter set with `KPG_PYTHON`.

`code/code1_IQR_calculation.py` is the interactive single-record processor. It requires an Excel file with `ID`, `Age` and `Index`, one operational record per file, and matching units for ages and bin edges. Use the compiled-inventory entry above for the supplied full-database calculation.

## Plot saved results

These commands read supplied summaries and do not rerun the uncertainty simulations or rarefaction.

```sh
Rscript code/modules/Code02_smoothing_cluster_bootstrap/plots/plot_Fig2A_background_comparison_Arial_20261009.R
Rscript code/modules/Code02_smoothing_cluster_bootstrap/plots/plot_FigS3_final_from_locked_inputs.R
Rscript code/modules/Code02_smoothing_cluster_bootstrap/plots/plot_ancillary_FigS10_from_locked_inputs.R
Rscript code/modules/Code02_smoothing_cluster_bootstrap/plots/plot_EV_proxy_correlations_final.R
Rscript code/modules/Code03_sampling_density_rarefaction/plots/plot_FigS4_Spearman_left_aligned_20261008.R
Rscript code/modules/Code03_sampling_density_rarefaction/plots/plot_FigS5_final.R
Rscript code/modules/Code04_measurement_uncertainty/02_code_revised/plot_FigS6_final.R
Rscript code/modules/Code05_chronology_uncertainty/02_code_revised/plot_FigS7_match_user_FigS6_Arial_20261009.R
Rscript code/modules/Code06_joint_uncertainty/02_code_revised/plot_FigS8_joint_conservative_measurement_chronology_uncertainty.R
Rscript code/modules/Code08_EV_proxy_biota_Pearson/Plot_Fig3.R
Rscript code/modules/Code09_mean_state_biota_Pearson/Plot_FigS13.R
Rscript code/modules/Code10_ADF_Johansen_VECM/Plot_Fig4.R
python3 code/modules/Code07_latitude/Code07_latitude_analysis.py --plot-only
```

New outputs use `generated/`, component `outputs/`, or the specified output directories. Plot-only commands read the supplied saved summaries rather than automatically substituting newly calculated summaries.

## Inputs and scope

See [INPUTS.md](INPUTS.md). Biological rates are supplied inputs, not re-estimated with PyRate here. Nannofossil rates use a supplied 100-kyr grid derived from approximately 1-Myr source data.

The helpers in `code/auxiliary/` smooth supplied SST and nannofossil extinction-percentage series; they are not complete upstream calibration or biological-rate estimation pipelines. Fig. 5 uses a supplied image for biological/event artwork, not a reconstruction from fossil occurrences.

The repository contains code, inputs and required reference summaries, not manuscript DOCX files. `FILE_MANIFEST.tsv` lists file checksums.
