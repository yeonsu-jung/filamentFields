# %%
import numpy as np
import matplotlib.pyplot as plt
import importlib
import filamentFields
import re

# %%
data_path = 'data/output.txt'

with open(data_path,'r') as f:
    lines = f.readlines()
    
# %%
# Parse the provided text file to extract the 3D arrays associated with each label.

file_path = 'data/output.txt'

# Read the file content
with open(file_path, 'r') as file:
    data = file.readlines()

# Initialize variables
arrays_per_label = {}
current_label = None
current_array = []

# Parse the file
for line in data:
    line = line.strip()
    if line.startswith('label'):
        # If there's an ongoing label, store its data
        if current_label is not None:
            arrays_per_label[current_label] = current_array
        # Reset for the new label
        current_label = int(line.split()[1])
        current_array = []
    else:
        # Append coordinates to the current array
        coords = list(map(int, line.split()))
        current_array.append(coords)

# Don't forget to save the last read label
if current_label is not None:
    arrays_per_label[current_label] = current_array

# %% filamentFields gets an input as a list of nx3 arrays.
list_of_arrays = [np.array(arrays_per_label[label]) for label in arrays_per_label]
fig,ax=plt.subplots(subplot_kw={'projection':'3d'})
for r in list_of_arrays:
    ax.plot(r[:,0],r[:,1],r[:,2])
plt.show()
    
# %%
fF = filamentFields.filamentFields(list_of_arrays)

R_omega = 30
# %%
local_edges_labels = fF.analyze_local_volume_from_precomputed(np.array([0,0,0]), R_omega, 0.01)

# %%
all_nodes = np.vstack(list_of_arrays)
xlim = [np.min(all_nodes[:,0]),np.max(all_nodes[:,0])]
ylim = [np.min(all_nodes[:,1]),np.max(all_nodes[:,1])]
zlim = [np.min(all_nodes[:,2]),np.max(all_nodes[:,2])]

num_1dgrid_points = 30
mg = np.meshgrid(np.linspace(xlim[0],xlim[1],num_1dgrid_points),
                    np.linspace(ylim[0],ylim[1],num_1dgrid_points),
                    np.linspace(zlim[0],zlim[1],num_1dgrid_points))

query_points = np.vstack([mg[0].flatten(),mg[1].flatten(),mg[2].flatten()]).T

# %%
import time

entanglement_fields_local = np.zeros(num_1dgrid_points**3)
start_time = time.time()
for i,query_point in enumerate(query_points):
    local_edges_labels = fF.analyze_local_volume(query_point, R_omega, 0.01)
    local_entanglement = fF.return_entanglement()
    entanglement_fields_local[i] = local_entanglement
print(f'Time taken: {time.time() - start_time}')

entanglement_fields_precomputed = np.zeros(num_1dgrid_points**3)
start_time = time.time()
fF.update_filament_nodes_list(list_of_arrays)
fF.precompute(R_omega)
for i,query_point in enumerate(query_points):
    local_edges_labels = fF.analyze_local_volume_from_precomputed(query_point, R_omega, 0.01)
    local_entanglement = fF.return_entanglement()
    entanglement_fields_precomputed[i] = local_entanglement
print(f'Time taken: {time.time() - start_time}')
    
# %%
entanglement_image_local = entanglement_fields_local.reshape(num_1dgrid_points,num_1dgrid_points,num_1dgrid_points)
entanglement_image_precomputed = entanglement_fields_precomputed.reshape(num_1dgrid_points,num_1dgrid_points,num_1dgrid_points)

print(f'Are they exactly same?: not yet...')
print(f'Fraction of identical elements: {np.sum(entanglement_image_local == entanglement_image_precomputed)/np.prod(entanglement_image_local.shape)}')
print(f'But they are very close to each other...')
print(f'Mean absolute difference: {np.mean(np.abs(entanglement_image_local - entanglement_image_precomputed))}')
print(f'Mean entanglement value: {np.mean(entanglement_image_local)}')

# %%
projected_along_z_local = np.sum(entanglement_image_local,axis=2)
projected_along_z_precomputed = np.sum(entanglement_image_precomputed,axis=2)

fig,ax=plt.subplots(1,2)
plt.colorbar(ax[0].imshow(projected_along_z_local))
plt.colorbar(ax[1].imshow(projected_along_z_precomputed))
plt.show()

# %%
projected_along_y_local = np.sum(entanglement_image_local,axis=1)
projected_along_y_precomputed = np.sum(entanglement_image_precomputed,axis=1)

fig,ax=plt.subplots(1,2)
plt.colorbar(ax[0].imshow(projected_along_y_local))
plt.colorbar(ax[1].imshow(projected_along_y_precomputed))
plt.show()



    

# %%
