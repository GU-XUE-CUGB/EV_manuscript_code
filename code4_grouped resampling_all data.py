import pandas as pd
import numpy as np
import os  # Used for checking if the directory exists

def resample_and_save_to_csv(input_file, output_csv_file, resample_size=3, resample_count=1000):
    """
    Perform resampling on tabular data, generate resampled points, and export the results to a CSV file.

    Parameters:
    - input_file: Input CSV file path containing 'Midpoint' and 'IQR' columns
    - output_csv_file: Output CSV file path to store the resampling results
    - resample_size: The sample size for each resampling
    - resample_count: The number of resampling iterations for each 'Midpoint' value
    """
    # Read the data
    data = pd.read_csv(input_file)
    x_column = 'Midpoint'  # Assuming the column for x is named 'Midpoint'
    y_column = 'IQR'  # Assuming the column for y is named 'IQR'

    # Store the resampling results
    resampled_x = []
    resampled_y = []

    # Perform resampling for each unique 'Midpoint'
    for x_value, group in data.groupby(x_column):
        if len(group) >= resample_size:
            for _ in range(resample_count):
                sampled_y = group[y_column].sample(n=resample_size, replace=True).mean()
                resampled_x.append(x_value)
                resampled_y.append(sampled_y)

    # Convert to NumPy arrays
    resampled_x = np.array(resampled_x)
    resampled_y = np.array(resampled_y)

    # Save the resampling results to the CSV file
    if output_csv_file:
        # Check if the directory exists
        output_dir = os.path.dirname(output_csv_file)
        if not os.path.exists(output_dir):
            print(f"Path {output_dir} does not exist, creating it...")
            os.makedirs(output_dir)

        # Save the resampled data as a DataFrame and export it to CSV
        result_df = pd.DataFrame({
            'Midpoint': resampled_x,
            'IQR': resampled_y
        })
        
        # Confirm file path and save the results
        print(f"Saving resampling results to: {output_csv_file}")
        result_df.to_csv(output_csv_file, index=False)
        print(f"Resampling results have been saved to {output_csv_file}")

# Example usage
input_file = ''  # Input CSV file path
output_csv_file = ''  # Output CSV file path

# Execute the function
resample_and_save_to_csv(input_file=input_file, output_csv_file=output_csv_file)
