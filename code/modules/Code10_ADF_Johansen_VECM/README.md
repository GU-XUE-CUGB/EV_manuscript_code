# Code 6 — ADF, Johansen, VECM and directional tests

The numerical entry points are in the root `code/` directory:

- `Rscript code/code6_Johansen_rate_screening.R`: three-rate screening; it also executes the two extinction models.
- `Rscript code/code6_ADF_Johansen_VECM.R`: the two extinction models alone.

Run from the repository root with `urca` installed. Ages decrease from older to younger. The extinction script embeds 48 elasmobranch pairs (50-kyr EV) and 24 bivalve pairs (100-kyr EV). Screening reads the current Code 5 pair table.

ADF uses `ur.df(type="none", selectlags="AIC", lags=1)`. Johansen uses `ca.jo(type="trace", ecdet="none", K=2, spec="longrun")`; rank is assessed sequentially at 5% and 10%. Eligible rank-1 systems use `cajorls(r=1)`. EV is first and its coefficient is normalized to +1.

ECT is extracted by equation and `ect1`. Short-run tests restrict the cross-variable lagged difference; joint tests restrict ECT plus that cross-lag. These tests measure predictive information, not physical causation.

New outputs go to `code/results/generated/code6_current_dataset/`. The saved statistics used for Fig. 4 are in `code/results/code6_current_dataset/`.

Plot command: `Rscript code/modules/Code10_ADF_Johansen_VECM/Plot_Fig4.R`. Fig. 4A uses the three-rate trace statistics; Fig. 4B uses the two extinction directional tests. The source uses Helvetica.
