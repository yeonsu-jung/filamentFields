import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data from CSV
file_path = '/Users/romyaran/volume_fraction_results_z100.csv'  # Replace with your CSV file path
data = pd.read_csv(file_path)

# Ensure the CSV file has columns 'x', 'y', 'z', 'vf'
if set(['x', 'y', 'z', 'vf']).issubset(data.columns):
    x = data['x']
    y = data['y']
    z = data['z']
    vf = data['vf']
    
    # Create 3D scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x, y, z, c=vf, cmap='magma', marker='o')
    
    # Add color bar which maps values to colors
    cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
    cbar.set_label('vf')
    
    # Set labels
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title('3D Scatter Plot with Volume Fraction')
    
    # Show plot
    plt.show()
else:
    print("CSV file must contain 'x', 'y', 'z', and 'vf' columns.")

import pandas as pd
import numpy as np

# Load the data from CSV files
file1 = '/Users/romyaran/volume_fraction_results_z100.csv'  # Replace with the path to your first CSV file
file2 = '/Users/romyaran/Desktop/Summer_2024/yeonsu_neuron_linkage_code_cpp/filamentFields-main/cross_sections/z/100/entanglement_results_z100.csv'  # Replace with the path to your second CSV file

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Merge the two dataframes on x, y, z coordinates
merged_df = pd.merge(df1, df2, on=['x', 'y', 'z'])

# Extract the fields of interest
field1 = merged_df['entanglement'].values
field2 = merged_df['vf'].values

# Calculate the mean of each field
mean_field1 = np.mean(field1)
mean_field2 = np.mean(field2)

# Calculate the covariance between the two fields
covariance = np.mean((field1 - mean_field1) * (field2 - mean_field2))

# Calculate the standard deviation of each field
std_field1 = np.std(field1)
std_field2 = np.std(field2)

# Calculate the correlation coefficient
correlation_coefficient = covariance / (std_field1 * std_field2)

print("Correlation Coefficient:", correlation_coefficient)
