import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data from CSV
file_path = '/Users/romyaran/Desktop/Summer_2024/yeonsu_neuron_linkage_code_cpp/filamentFields-main/cross_sections/z/100/entanglement_results_z100.csv'  # Replace with your CSV file path
data = pd.read_csv(file_path)

# Ensure the CSV file has columns 'x', 'y', 'z', 'entanglement'
if set(['x', 'y', 'z', 'entanglement']).issubset(data.columns):
    x = data['x']
    y = data['y']
    z = data['z']
    entanglement = data['entanglement']
    
    # Create 3D scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x, y, z, c=entanglement, cmap='viridis', marker='o')
    
    # Add color bar which maps values to colors
    cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
    cbar.set_label('Entanglement')
    
    # Set labels
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title('3D Scatter Plot with Entanglement')
    
    # Show plot
    plt.show()
else:
    print("CSV file must contain 'x', 'y', 'z', and 'entanglement' columns.")
