# Code 3 — sampling density and rarefaction

Run `python3 run.py` for 10,000 rarefaction iterations. Required packages: NumPy and pandas. Outputs are written under `outputs/formal_10000/`; paths resolve from this folder.

The analysis uses 49 standard 50-kyr bins (67.20–64.80 Ma) and the bin assignments supplied in the retained-measurement table. Eligible observed record-bins contain at least 2 measurements and have positive IQR.

Fixed-record rarefaction samples 3, 5 or 7 records without replacement per bin, with seeds 20260824, 20260826 and 20260828. Joint rarefaction samples 5 records per bin and 3, 4 or 5 measurements per record-bin without replacement. Seed 20260821 is reset for each measurement-depth setting. Record-wide z-scores are retained. Zero IQRs arising within resampled subsets remain in the analysis.

Sampling-density associations use original bin-level series and older-to-younger first differences. Paired circular moving-block bootstrap intervals use block lengths 3, 4 and 5 (primary: 4), 10,000 resamples and linearly interpolated percentiles. Seeds are defined in the source.

Run `Rscript run_density_nominal_P.R` for the eight two-sided nominal Spearman P values using `cor.test(..., exact=FALSE)`.

Plot commands:

- `Rscript plots/plot_FigS4_Spearman_left_aligned_20261008.R`
- `Rscript plots/plot_FigS5_final.R`

Plot scripts use the supplied summaries in `plot_snapshot_inputs/` and write under `outputs/plots/`. Fig. S5 display curves use `npreg::ss(spar=0.2)`. R plotting packages include readr, dplyr, tidyr, ggplot2, npreg, gridExtra, svglite and ragg; the font is Arial.
