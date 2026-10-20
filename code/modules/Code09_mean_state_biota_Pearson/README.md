# Code 5 — mean environmental states and biological Pearson correlations

Run `Rscript Supplementary_Code09.R`; draw Fig. S13 with `Rscript Plot_FigS13.R`. Analysis uses base R; plotting requires grid, ragg and svglite. Outputs go to `results/`.

Inputs contain 10,881 retained raw measurements in 70 operational records. Each record is standardized using its mean and sample SD. Within-bin mean state is the arithmetic mean of individual standardized measurements, not a mean of record IQRs or record means. Global state includes retained measurements in the pooled pCO₂ record.

Fifty-kyr bins start at 64.775 Ma; 100-kyr bins start at 64.75 Ma, using the age conditions implemented in the source. No IQR eligibility mask is applied to mean-state measurements.

Environmental and biological ages are matched exactly where available. Otherwise `stats::approx(rule=1)` is restricted to two finite neighboring bins separated by no more than one bin width. There is no endpoint carry or interpolation across missing environmental bins. Biological grids and rate values are supplied inputs.

The Fig. S13 table contains 45 two-sided Pearson tests for Global state and four proxy states. Nine additional pCO₂ tests are output separately, giving 54 rows in the complete table. Rate standardization uses `scale`. P values are nominal and n is the finite paired count. The supplied plot source uses Helvetica.
