# Code 1 — EV reconstruction

Run `python3 run.py`. Required packages: NumPy and pandas. Paths resolve from this folder; outputs are written to `outputs/generated/`.

Input: the compiled measurement inventory, including raw values and record information. The script applies the specified 64-value Wang-derived pCO₂ exclusion before standardization. The pooled pCO₂ operational record contains 67 retained measurements.

Each operational record is standardized using its record-wide mean and sample SD (`ddof=1`). Measurements are assigned to the configured 50-, 100-, 150- and 200-kyr bins. Within-record IQR uses linearly interpolated quartiles and requires at least 2 measurements per record-bin. Positive record-bin IQRs contribute to the record-weighted arithmetic mean EV.

The background value is the equal-weight mean of the pre-window mean (66.95–67.20 Ma, 6 bins) and post-window mean (64.80–65.45 Ma, 14 bins): `0.34503257945945265`.

`lib/` contains the calculation functions. `reference_results/` contains saved comparison tables. Source-workbook extraction into the compiled inventory is a required upstream input step. The script generates numerical reconstruction tables, not the complete manuscript artwork.
