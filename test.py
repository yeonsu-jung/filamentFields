# %%
import numpy as np
import matplotlib.pyplot as plt
import importlib
import filamentFields
importlib.reload(filamentFields)

filaments = []
for i in range(1500):
    filaments.append(np.random.randn(10,3))

R_omega = 1
num_repeat = 1
    
fF = filamentFields.filamentFields(filaments)
fF.precompute(R_omega)
# %%
fig,ax=plt.subplots()
ax.plot([],[])

# %%
fF.compute_total_linking_matrix()
total_lk = fF.return_total_linking_matrix()
np.sum(np.abs(total_lk[~np.isnan(total_lk)]))
# %%

import time
start = time.time()
# fF.compute_total_linking_matrix() # to be fair
for i in range(num_repeat):
    fF.analyze_local_volume_from_precomputed(np.array([0,0,0]), R_omega,0.01)
print(f'Time taken: {time.time() - start}')
print(fF.return_entanglement())
print(fF.return_number_of_labels())
# print(fF.return_orientational_order_parameter())
# %%
start = time.time()
for i in range(num_repeat):
    fF.analyze_local_volume(np.array([0,0,0]), R_omega,0.01)
print(f'Time taken: {time.time() - start}')
print(fF.return_entanglement())
# print(fF.return_number_of_labels())
# print(fF.return_orientational_order_parameter())
# %%





# %%
filaments = []
for i in range(1000):
    x = np.cumsum(np.random.randn(100))
    y = np.cumsum(np.random.randn(100))
    z = np.cumsum(np.random.randn(100))
    x = np.convolve(x, np.ones(5)/5, mode='valid')
    y = np.convolve(y, np.ones(5)/5, mode='valid')
    z = np.convolve(z, np.ones(5)/5, mode='valid')        
    filaments.append(np.vstack([x,y,z]).T)
    

# %%
fig,ax=plt.subplots(subplot_kw={'projection':'3d'})
for r in filaments:
    ax.plot(r[:,0],r[:,1],r[:,2])
    
# %%
filaments = np.array(filaments)
xlim = [np.min(filaments[:,:,0]),np.max(filaments[:,:,0])]
ylim = [np.min(filaments[:,:,1]),np.max(filaments[:,:,1])]
zlim = [np.min(filaments[:,:,2]),np.max(filaments[:,:,2])]

num_1dgrid_points = 10
mg = np.meshgrid(np.linspace(xlim[0],xlim[1],num_1dgrid_points),np.linspace(ylim[0],ylim[1],num_1dgrid_points),np.linspace(zlim[0],zlim[1],num_1dgrid_points))

query_points = np.vstack([mg[0].flatten(),mg[1].flatten(),mg[2].flatten()]).T
# %%



# %%

entanglement_results = np.zeros(num_1dgrid_points**3)
number_of_labels = np.zeros(num_1dgrid_points**3)
for i in range(num_1dgrid_points**3):
    local_edges = fF.analyze_local_volume(query_points[i], 5,0.01)
    
    entanglement_results[i] = fF.return_entanglement()
    number_of_labels[i] = fF.return_number_of_labels()
    
    
    
# %%    
fig,ax=plt.subplots(subplot_kw={'projection':'3d'})
for i in range(len(local_edges)):
    edge = local_edges[i]
    p1 = edge[:3]
    p2 = edge[3:]
    ax.plot([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]])

    
# %%
fig,ax=plt.subplots(subplot_kw={'projection':'3d'})
ax.scatter(query_points[:,0],query_points[:,1],query_points[:,2],c=entanglement_results)

# %%
plt.hist(entanglement_results[entanglement_results>0],bins=100)

# %%
plt.hist(number_of_labels[number_of_labels>0],bins=100)
# %%
n_image = number_of_labels.reshape(num_1dgrid_points,num_1dgrid_points,num_1dgrid_points)
projected_image = np.sum(n_image,axis=0)

fig,ax=plt.subplots()
ax.imshow(projected_image)

# %%
fF = filamentFields.filamentFields(filaments)
all_nodes = fF.return_all_nodes()

# %%
query_points = np.array([0,0,0])
R_omega = 2
local_edges = fF.analyze_local_volume(query_points, R_omega, 0.01)
# %%
fig,ax=plt.subplots(subplot_kw={'projection':'3d'})
for i in range(len(local_edges)):
    edge = local_edges[i]
    p1 = edge[:3]
    p2 = edge[3:]
    ax.plot([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]])
# %%
# my sampling
filament_edges_list = fF.return_filament_nodes_list()
all_edges = fF.return_all_edges()
# %%
len(filament_edges_list)
fig,ax=plt.subplots(subplot_kw={'projection':'3d'})
for edge in all_edges:
    p1 = edge[:3]
    p2 = edge[3:]
    
    if (np.linalg.norm(p1 - query_points) < R_omega ) & (np.linalg.norm(p2 - query_points) < R_omega):        
        ax.plot([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]])
    

# %%

# %%
num_rods = 100
def linspace3(start, stop, num):
    return np.vstack([np.linspace(start[i], stop[i], num) for i in range(3)]).T 

single_edge = linspace3([0,0,0],[0,0,1],10)

mg = np.meshgrid(np.linspace(0,1,num_rods),np.linspace(0,1,num_rods))
mg = [mg[0].flatten(),mg[1].flatten()]
aligned_filaments = []
for i in range(num_rods**2):
    translator = np.array([mg[0][i],mg[1][i],0]).T
    aligned_filaments.append(single_edge + translator)
    
# fig,ax=plt.subplots(subplot_kw={'projection':'3d'})
# for i in range(100):
#     ax.plot(aligned_filaments[i][:,0],aligned_filaments[i][:,1],aligned_filaments[i][:,2])
    
# %%
fF = filamentFields.filamentFields(aligned_filaments)
# %%
fF.analyze_local_volume(np.array([0,0,0]), 1,0.01)
