# Code 5 — environmental variability and biological Pearson correlations

Run `Rscript Supplementary_Code08.R --current-inputs`; draw Fig. 3 with `Rscript Plot_Fig3.R`. Analysis uses base R; plotting additionally requires ggplot2, patchwork, ragg, svglite and scales. Numerical outputs go to `results/generated/current_dataset/`; the plot reads supplied summaries in `results/current_dataset/`.

The supplied age-aligned environmental–biological pair table and environmental bin table are the inputs. Elasmobranchs use 50-kyr EV; bivalves and calcareous nannofossils use 100-kyr EV. The supplied nannofossil grid is derived from approximately 1-Myr source data.

The script retains supplied EV values aligned between neighboring observed bins, excludes endpoint carry and spans longer than one bin width, and creates no additional interpolation.

The main output contains 45 full-transition Pearson tests: five environmental metrics × three taxa × three biological variables. Two-sided nominal P values use `stats::cor.test(method="pearson")`; n is the finite paired count.

`Rscript Ancillary_stage_Pearson.R` produces the supplied stage-analysis scope: >66.0 Ma, 66.0–65.5 Ma, <65.5 Ma and full transition. Outputs go to `results/ancillary_existing_stages/`.

Fig. 3 uses three circle-size classes (P≥0.05, 0.005≤P<0.05, P<0.005), no r text and a symmetric −0.7 to +0.7 color scale. The supplied plotting source uses Helvetica. The additional pCO₂ column uses `Supplementary_pCO2_variability.R --current-inputs`.
