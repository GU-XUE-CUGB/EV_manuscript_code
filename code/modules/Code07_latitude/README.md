# Code 7 — latitude-specific EV and site-balanced rarefaction

Run `python Code07_latitude_analysis.py --no-plots` for numerical outputs; omit `--no-plots` to generate figures as well. After numerical outputs are available, `--plot-only` draws figures without resampling.

Required Python packages: NumPy and openpyxl. R plotting packages: ggplot2, npreg, patchwork, jsonlite, svglite, png, ragg and grid. Inputs and outputs are relative to this folder; results go to `output/`.

Descriptive latitude-specific EV uses positive 50-kyr record-bin IQRs and record weighting. The pooled pCO₂ record has no unique geographic site and is excluded. Latitude and canonical-site assignments come from the supplied mapping tables.

For site-balanced comparisons, records within a site and bin are first averaged. Both latitude groups must contain at least 2 independent sites. Each group supplies k=min(n_low,n_mid_high) sites without replacement. There are 10,000 iterations. Seed 20260929 uses NumPy PCG64 generators with SeedSequence.spawn: one baseline stream and one per sorted canonical site. Phase summaries and leave-one-site-out comparisons use the rules in the source.

Both descriptive EV curves use `npreg::ss(spar=0.25)`. Site-balanced differences and interval endpoints are not smoothed. K–Pg is 66.016 Ma.

Fig. 5 biological and event blocks are read from the supplied raster asset; the environmental contrast is vector. The supplied plot layout appends that contrast as panel D. Figs. S14 and S15 are drawn from numerical outputs.
