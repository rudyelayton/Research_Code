# -*- coding: utf-8 -*-
"""
Created on Mon May 11 18:38:22 2020

@author: rudye layton
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Input 
ndata = 200
pi = np.pi
sq3 = np.sqrt(3)
a = 1
s = 0
eps = 0

# Function
def w(kx,ky):
    W = np.sqrt(1.0 + 4.0*np.cos(0.5*sq3*kx*a)*np.cos(0.5*ky*a) + 4*(np.cos(0.5*ky*a))**2)
    return W

def M(kx,ky):
    M = ((np.cos(0.5*ky*a)*np.cos(0.5*kx*sq3*a) - np.cos(ky*a))**2) + (3*(np.sin(0.5*ky*a))**2*(np.sin(0.5*sq3*kx*a))**2)
    return M

# Brillouin Zone
kbzx1 = np.linspace(-2*pi*sq3/3, -2*pi*sq3/3, ndata)
kbzy1 = np.linspace(-4*pi*sq3/9, 4*pi*sq3/9, ndata)

kbzx2 = np.linspace(2*pi*sq3/3, 2*pi*sq3/3, ndata)
kbzy2 = np.linspace(-4*pi*sq3/9, 4*pi*sq3/9, ndata)

kbzx3 = np.linspace(2*pi*sq3/3, 0, ndata)
kbzy3 = np.linspace(-4*pi*sq3/9, -8*pi*sq3/9, ndata)

kbzx4 = np.linspace(-2*pi*sq3/3, 0, ndata)
kbzy4 = np.linspace(-4*pi*sq3/9, -8*pi*sq3/9, ndata)

kbzx5 = np.linspace(-2*pi*sq3/3, 0, ndata)
kbzy5 = np.linspace(4*pi*sq3/9, 8*pi*sq3/9, ndata)

kbzx6 = np.linspace(2*pi*sq3/3, 0, ndata)
kbzy6 = np.linspace(4*pi*sq3/9, 8*pi*sq3/9, ndata)

bzx = np.concatenate((kbzx2,kbzx3,kbzx4,kbzx1,kbzx5,kbzx6))
bzy = np.concatenate((kbzy2,kbzy3,kbzy4,kbzy1,kbzy5,kbzy6))

# Generate point for 3D point
kmesh = np.linspace(-2*pi, 2*pi, ndata)
kx , ky = np.meshgrid(kmesh, kmesh)

MT = (M(kx,ky)/w(kx,ky))**2

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Matrix Transition', fontsize=20)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.set_xlabel('kx', fontsize=13)
ax.set_ylabel('ky', fontsize=13)
ax.set_zlabel('M^2', fontsize=13)
Ec = ax.plot_surface(kx, ky, MT, cmap='viridis')
bz = ax.plot(bzx, bzy, 20, color='red')
plt.show()
print('Saving matrix transition probability as photo')
plt.savefig('GrapheneMatrixTransition.png')