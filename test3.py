# %%
import numpy as np
import os
from pathlib import Path

pth = Path('/scratch/08871/tg881707/fil_data/filaments_2048')

def compute_filament_length(filament):
    return np.sum(np.sqrt(np.sum(np.diff(filament,axis=0)**2,axis=1)))

def parse_file_to_numpy_arrays(file_path):
    arrays = []
    current_array = []
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('label'):
                if current_array:
                    arrays.append(np.array(current_array, dtype=int))
                    current_array = []
            else:
                current_array.append([int(x) for x in line.split()])
                
    if current_array:
        arrays.append(np.array(current_array, dtype=int))
        
    return arrays


filaments = []
iter = 0
for file_path in list(pth.glob('filament*'))[30:32]:
    filaments.extend(parse_file_to_numpy_arrays(file_path))
    iter += 1
    # to read a small number of filament; to be removed later    
    if iter > 1:
        break

filament_length_list = []
for filament in filaments:
    filament_length_list.append(compute_filament_length(filament))
avg_filament_length = np.mean(filament_length_list)
    
print(f'Number of filaments: {len(filaments)}')
print(f'Average length of filament: {avg_filament_length}')

from matplotlib import pyplot as plt
fig,ax=plt.subplots(1,1,subplot_kw={'projection':'3d'})
for filament in filaments:
    ax.plot(filament[:,0],filament[:,1],filament[:,2])
plt.savefig('test3.png')

import filamentFields
fF = filamentFields.filamentFields(filaments)

all_filaments = np.vstack(filaments)
xlim = np.min(all_filaments[:,0]),np.max(all_filaments[:,0])
ylim = np.min(all_filaments[:,1]),np.max(all_filaments[:,1])
zlim = np.min(all_filaments[:,2]),np.max(all_filaments[:,2])
print((xlim,ylim,zlim))

R_omega = 32
num_grid = 30
mg = np.meshgrid(np.linspace(xlim[0],xlim[1],num_grid),np.linspace(ylim[0],ylim[1],num_grid),np.linspace(zlim[0],zlim[1],num_grid))
sampling_points = np.array([mg[0].flatten(),mg[1].flatten(),mg[2].flatten()]).T


import time
start_time = time.time()
filament_fields_result = fF.analyze_local_volume_over_domain(sampling_points, R_omega, 0.01)    
print(f'Elapsed: {time.time() - start_time} sec')

# n_field = filament_fields_result[:,0]
e_field = filament_fields_result[:,3]
max_e = np.max(e_field[~np.isnan(e_field)])
print(f'Max entanglement: {max_e}')

e_volume = e_field.reshape(num_grid,num_grid,num_grid).copy()
e_volume[np.isnan(e_volume)] = 0
e_proj_z = np.mean(e_volume,axis=0)

plt.close('all')
plt.imshow(e_proj_z)
plt.axis('equal')
plt.savefig('e_proj_z.png')