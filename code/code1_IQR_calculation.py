import pandas as pd
import numpy as np

# Input: one final retained operational record per file; ages and bin edges use the same units.

# Step 1: Read EXCEL data and perform preliminary statistical calculations

# 1.1 Read data
file_path = input("Enter the EXCEL file path: ")
data = pd.read_excel(file_path)

# 1.2 Validate Age data
if not {'ID', 'Age', 'Index'}.issubset(data.columns):
    raise ValueError("Input requires ID, Age and Index columns")
if not np.isfinite(data[['Age', 'Index']].to_numpy(dtype=float)).all():
    raise ValueError("Age and Index must contain finite numeric values")
for record_column in ['record_id', 'Record']:
    if record_column in data and data[record_column].nunique() > 1:
        raise ValueError("Use one operational record per file, or the whole-database reconstruction module")

# 1.3 Calculate Z-Score
M = data['Index'].mean()
SD = data['Index'].std(ddof=1)
if not np.isfinite(SD) or SD <= 0:
    raise ValueError("Record-wide sample SD must be positive")
data['Z-Score'] = (data['Index'] - M) / SD

# Step 2: Define and assign TimeBins

def assign_timebin(age):
    for i, (start, end) in enumerate(bins):
        if start <= age < end:
            return f"TimeBin{i+1}"
    return None

# 2.1 Get user input for time range and bin size
Start = float(input("Enter Start Age: "))
End = float(input("Enter End Age: "))
Bin = float(input("Enter TimeBin size: "))
if not np.isfinite([Start, End, Bin]).all() or Start <= End or Bin <= 0:
    print("Invalid Age Range")
    exit()

# 2.2 Define time intervals and assign bins
bins = []
midpoints = []
n = round((Start - End) / Bin)
if not np.isclose((Start - End) / Bin, n, rtol=0, atol=1e-8):
    raise ValueError("Age range must contain an integer number of fixed-width bins")
for i in range(n):
    bin_start, bin_end = End + i * Bin, End + (i + 1) * Bin
    bins.append((bin_start, bin_end))
    midpoints.append((bin_start + bin_end) / 2)
data['TimeBin'] = data['Age'].apply(assign_timebin)

# 2.3 Count the number of data points in each TimeBin
timebin_counts = []
for bin_start, bin_end in bins:
    count = len(data[(data['Age'] >= bin_start) & (data['Age'] < bin_end)])
    timebin_counts.append(count)

# Step 3: Calculate IQR(INC)

# 3.1 Compute IQR(INC)
iqr_values_inc = []
for bin_start, bin_end in bins:
    subset = data[(data['Age'] >= bin_start) & (data['Age'] < bin_end)]['Z-Score'].sort_values()
    if len(subset) >= 2:
        Q1_inc = subset.quantile(0.25, interpolation='linear')
        Q3_inc = subset.quantile(0.75, interpolation='linear')
        iqr_values_inc.append(Q3_inc - Q1_inc)
    else:
        iqr_values_inc.append(None)

# 3.2 Merge IQR data
iqr_df = pd.DataFrame({'TimeBinIndex': [f"TimeBin{i+1}" for i in range(len(bins))],
                       'Time midpoint for plotting': midpoints,
                       'IQR(INC)': iqr_values_inc})
data = pd.merge(data, iqr_df, left_on='TimeBin', right_on='TimeBinIndex', how='left')

# Step 4: Output data to EXCEL file

# Save data to an EXCEL file
output_path = input("Enter output EXCEL file path: ")

with pd.ExcelWriter(output_path) as writer:
    
    # Create TimeBin_Information DataFrame
    timebin_info = pd.DataFrame({
        'TimeBin': [f"TimeBin{i + 1}" for i in range(len(bins))],
        'Bin_end': [bin_end for _, bin_end in bins],
        'Bin_start': [bin_start for bin_start, _ in bins],
        'Bin_midpoint': midpoints,
        'Data Counts': timebin_counts,
        'IQR(INC)': iqr_values_inc
    })

    # Save TimeBin_Information to an EXCEL sheet
    timebin_info.to_excel(writer, sheet_name='TimeBin_Information', index=False)

    # Save processed data to "Index_Calculation" sheet
    data[['ID', 'Age', 'Index', 'Z-Score', 'TimeBin', 'Time midpoint for plotting', 'IQR(INC)']].to_excel(
        writer, sheet_name='Index_Calculation', index=False)

print("Data has been saved to", output_path)

# Whole-database reconstruction: python modules/Code01_EV_reconstruction/run.py
# Zero IQRs remain in the output; the main 50-kyr EV mean excludes them.
