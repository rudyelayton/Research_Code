# -*- coding: utf-8 -*-
"""
Created on Mon May 11 12:10:54 2020

@author: ASUS
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Input 
ndata = 200
pi = np.pi
a = 1
t = -1
d = 0.5

# Generate point
kmesh = np.linspace(-1, 1, ndata)
kx , ky = np.meshgrid(kmesh, kmesh)
k = np.sqrt(kx**2 + ky**2)

# Function
miulit = d/np.sqrt(d**2 + 3*(t**2)*(a**2)*(k**2))

# Output 3D Polarization
write = np.column_stack((kx.flatten(), miu.flatten(), miulit.flatten()))
print('Saving polarization data as excel')
np.savetxt('hBNPolarization.csv', write, delimiter = ",", header=("k, miu, miulit"))

# Plot Pseudospin Polarization
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Pseudospin Polarization $\Delta$ = 0.5', fontsize=20)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel('kx', fontsize=13)
ax.set_ylabel('ky', fontsize=13)
ax.set_zlabel('Polarization', fontsize=13)
Ec = ax.plot_surface(kx, ky, miulit, cmap='coolwarm')
cb = plt.colorbar(Ec)
plt.show()