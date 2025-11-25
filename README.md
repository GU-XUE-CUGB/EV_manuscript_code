# EV_manuscript_code
This repository contains all scripts used in the analyses for our study reconstructing environmental variability (EV) across the Cretaceous–Paleogene transition and evaluating its ecological consequences.
The code supports the full workflow of EV reconstruction, statistical analyses, robustness assessments, and comparisons with marine biotic turnover rates.

Repository Structure
1. code1_IQR_calculation
Calculates the interquartile range (IQR) for each paleoclimate proxy record.
2. code2_smoothing_spline
Applies a smoothing spline to the reconstructed EV time series. This smoothing is used exclusively for visualization, to clarify long-term EV trends.
3. code3_bootstrap_resampling_all_data
Performs bootstrap resampling on the 50-kyr time-binned global EV reconstruction.
This evaluates the robustness of the global EV pattern by repeatedly resampling EV values across bins.
4. code4_grouped_resampling_all_data
Conducts proxy-group bootstrap resampling on the 50-kyr time-binned global EV reconstruction.
5. code5_Pearson_analysis
Computes Pearson correlations for three categories of relationships:
Global EV vs. proxy-specific EV.
Global EV vs. biotic turnover rates of marine groups (elasmobranchs, bivalves, calcareous microfossils). 
Proxy-specific EV vs. biotic turnover rates of marine groups (elasmobranchs, bivalves, calcareous microfossils).
These analyses quantify the strength and significance of links between environmental variability and marine ecosystem disruption.
6. code6_ADF+Johansen_analysis+VECM
Performs stationarity testing (ADF), Johansen cointegration analysis, and Vector Error Correction Modeling (VECM) to evaluate the long-term equilibrium relationships between global EV and marine biotic turnover rates (elasmobranchs, bivalves).
This identifies whether EV can be interpreted as a long-term driver of extinction dynamics.
7. code7_bootstrap_resampling_mid-high_latitude
Bootstraps the IQR datasets from mid- to high-latitude sites, repeatedly drawing subsets of IQRs and reconstructing regional EV curves.
Used to evaluate the robustness of mid–high-latitude EV and its sensitivity to sample number.
8. code8_bootstrap_resampling_low_latitude
Bootstraps the IQR datasets from low-latitude sites, repeatedly resampling subsets of proxy IQRs to rebuild regional EV curves.
Used to assess the robustness of low-latitude EV patterns, facilitating comparison with mid–high-latitude regions.

Data Availability
The paleoclimate proxy datasets used for EV reconstruction are provided in the manuscript’s Supplementary Materials (Tables S2–S6).
These include global records of temperature, δ¹⁸O, δ¹³C, atmospheric pCO₂, and weathering indicators.
Biotic turnover rates (extinction, speciation, and net diversification) for elasmobranchs and bivalves were obtained from:
Guillaume Guinot & Fabien L. Condamine (2023), Science
DOI: 10.1126/science.abn2080
Calcareous nannoplankton data were obtained from:
Lowery, Bown, Fraass & Hull (2020), Annual Review of Earth and Planetary Sciences
DOI: 10.1146/annurev-earth-081619-052818
The interpolated nannoplankton dataset used in our analyses is provided in the manuscript (Table S8).

Purpose
This repository enables full reproducibility of:
EV reconstruction
Proxy-specific and regional variability analysis
Robustness assessments (global and latitudinal)
Correlation and cointegration analyses
Ecological comparisons with marine turnover rates

Users can directly follow the workflow or adapt individual components for related paleoclimate or macroevolutionary studies.
