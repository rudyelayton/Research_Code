# -*- coding: utf-8 -*-
"""
Created on Fri May  8 22:24:27 2020

@author: ASUS
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# generate dummy data (graphene tight binding band structure)
kvec = np.linspace(-np.pi,np.pi,101)
kx,ky = np.meshgrid(kvec,kvec)
E = np.sqrt(1+4*np.cos(3*kx/2)*np.cos(np.sqrt(3)/2*ky) + 4*np.cos(np.sqrt(3)/2*ky)**2)

# plot full dataset
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(kx,ky,E,cmap='viridis',vmin=-E.max(),vmax=E.max(),rstride=1,cstride=1)
ax.plot_surface(kx,ky,-E,cmap='viridis',vmin=-E.max(),vmax=E.max(),rstride=1,cstride=1)



# focus on Dirac cones
Elim = 1  #threshold
E[E>Elim] = np.nan

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(kx,ky,E,cmap='viridis',vmin=-Elim,vmax=Elim)
#ax.plot_surface(kx,ky,-E,cmap='viridis',vmin=-Elim,vmax=Elim)
ax.plot_surface(kx,ky,E,cmap='viridis',rstride=1,cstride=1,vmin=-Elim,vmax=Elim)
ax.plot_surface(kx,ky,-E,cmap='viridis',rstride=1,cstride=1,vmin=-Elim,vmax=Elim)

plt.show()