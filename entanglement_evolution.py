# %%
import numpy as np


pth = '/Users/yeonsu/GitHub/filamentFields/movie_frames/entanglement_over_time.csv'
data = np.loadtxt(pth, delimiter=',',skiprows=1)

time_steps = data[:,0]
entanglements = data[:,1]

import matplotlib.pyplot as plt
plt.figure()
plt.plot(time_steps, entanglements,'.-')
plt.xlabel('Time step')
plt.ylabel('Entanglement')
plt.title('Entanglement evolution over time')
plt.savefig('entanglement_evolution.png')
print('Saved entanglement_evolution.png')