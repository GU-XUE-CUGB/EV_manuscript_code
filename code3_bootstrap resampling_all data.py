import pandas as pd
import numpy as np

# Set input and output paths
input_file = ''  # Input CSV file path
output_csv_file = ''  # Output CSV file path

# Read the input CSV file
data = pd.read_csv(input_file)

# Assume the CSV file has 'Midpoint' and 'IQR' columns
Midpoint = data['Midpoint'].values
IQR = data['IQR'].values

# Set the number of resamples and the sample size for each resample
n_resamples = 1000
sample_size = 10

# Create empty lists to store all resampled points
resampled_midpoint_all = []
resampled_iqr_all = []

# Perform resampling
for _ in range(n_resamples):
    # Randomly select indices with replacement
    indices = np.random.choice(range(len(Midpoint)), size=sample_size, replace=True)
    resampled_midpoint = Midpoint[indices]
    resampled_iqr = IQR[indices]
    resampled_midpoint_all.extend(resampled_midpoint)
    resampled_iqr_all.extend(resampled_iqr)

# Save all resampled points into a DataFrame
resampled_df = pd.DataFrame({'resample_midpoint': resampled_midpoint_all, 'resample_iqr': resampled_iqr_all})

# Save the resampling results to the CSV file
resampled_df.to_csv(output_csv_file, index=False)

print(f"Resampling results have been saved to {output_csv_file}")
