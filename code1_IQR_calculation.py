import pandas as pd

# Step 1: Read EXCEL data and perform preliminary statistical calculations

# 1.1 Read data
file_path = input("Enter the EXCEL file path: ")
data = pd.read_excel(file_path)

# 1.2 Validate Age data
invalid_ids = [data['ID'].iloc[i+1] for i in range(len(data['Age'])-1) if data['Age'].iloc[i+1] < data['Age'].iloc[i]]
if invalid_ids:
    print("Invalid Age Information for IDs:", ', '.join(invalid_ids))
    exit()

# 1.3 Calculate Z-Score
M = data['Index'].mean()
SD = data['Index'].std()
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
if Start <= End:
    print("Invalid Age Range")
    exit()

# 2.2 Define time intervals and assign bins
bins = []
midpoints = []
n = 0
while End + n * Bin < Start:
    bins.append((End + n * Bin, End + (n + 1) * Bin))
    midpoints.append((End + n * Bin + End + (n + 1) * Bin) / 2)
    n += 1
if bins[-1][1] < Start:
    bins.append((bins[-1][1], Start))
    midpoints.append((bins[-1][1] + Start) / 2)
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
        Q1_inc = subset.quantile(0.25)
        Q3_inc = subset.quantile(0.75)
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
