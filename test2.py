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
local_edges_labels = fF.analyze_local_volume(np.array([0,0,0]), R_omega, 0.01)
# group by labels (7th column)
local_labels = np.unique(np.array(local_edges_labels)[:,6]).astype(int)
local_edges = []
for label in local_labels:
    edges = np.array(local_edges_labels)
    edges = edges[edges[:,6] == label,:6]
    local_edges.append(edges)
    
# %%
fig,ax=plt.subplots(subplot_kw={'projection':'3d'})
for i in range(len(local_edges)):
    edge = local_edges[i]
    node = np.vstack((edge[:,:3],edge[-1,3:]))
    ax.plot(node[:,0],node[:,1],node[:,2])
plt.show()
# %%
local_entanglement = fF.return_entanglement()
print(f'local entanglement: {local_entanglement}')

# %%
all_nodes = np.vstack(list_of_arrays)
xlim = [np.min(all_nodes[:,0]),np.max(all_nodes[:,0])]
ylim = [np.min(all_nodes[:,1]),np.max(all_nodes[:,1])]
zlim = [np.min(all_nodes[:,2]),np.max(all_nodes[:,2])]

num_1dgrid_points = 10
mg = np.meshgrid(np.linspace(xlim[0],xlim[1],num_1dgrid_points),
                    np.linspace(ylim[0],ylim[1],num_1dgrid_points),
                    np.linspace(zlim[0],zlim[1],num_1dgrid_points))

query_points = np.vstack([mg[0].flatten(),mg[1].flatten(),mg[2].flatten()]).T

# %%
entanglement_fields = np.zeros(num_1dgrid_points**3)
for i,query_point in enumerate(query_points):
    local_edges_labels = fF.analyze_local_volume(query_point, R_omega, 0.01)
    local_entanglement = fF.return_entanglement()
    entanglement_fields[i] = local_entanglement
    
# %%
entanglement_image = entanglement_fields.reshape(num_1dgrid_points,num_1dgrid_points,num_1dgrid_points)
projected_along_z = np.sum(entanglement_image,axis=2)
fig,ax=plt.subplots()
ax.imshow(projected_along_z)
plt.show()

# %%
projected_along_y = np.sum(entanglement_image,axis=1)
fig,ax=plt.subplots()
ax.imshow(projected_along_y)
plt.show()


    


