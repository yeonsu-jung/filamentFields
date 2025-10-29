# %%
import numpy as np
from matplotlib import pyplot as plt
import filamentFields

pth = '/Users/yeonsu/Downloads/nest_packing/positions/rod_coordinates_0000000000.csv'
data = np.loadtxt(pth, delimiter=',')

num_nodes = 9
all_rods = data.reshape(-1, num_nodes, 3)

# %%
fF = filamentFields.filamentFields(all_rods)

# %%
R_omega = 200.0
query_point = np.array([10.0, 10.0, 10.0])

# # quick non-precomputed analysis
local_edges = fF.analyze_local_volume(query_point, R_omega, 0.01)
print("Local edges shape:", local_edges.shape)
print("Entanglement:", fF.return_entanglement())
print("Number of labels:", fF.return_number_of_labels())

# %%


