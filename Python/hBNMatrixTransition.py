# -*- coding: utf-8 -*-
"""
Created on Mon May 11 18:38:22 2020

@author: ASUS
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
Kpr = -1
K = 1

# Function
def w(kx,ky):
    W = np.sqrt(1.0 + 4.0*np.cos(0.5*sq3*kx*a)*np.cos(0.5*ky*a) + 4*(np.cos(0.5*ky*a))**2)
    return W

def M1(kx,ky):
    M1 = (1 - (np.cos(1.5*kx*a/sq3)*(np.cos(0.5*ky*a) + sq3*Kpr*np.sin(0.5*ky*a))) + (np.cos(0.5*ky*a) + sq3*Kpr*np.sin(0.5*ky*a))**2)
    return M1

def M2(kx,ky):
    M2 = (1 - (np.cos(1.5*kx*a/sq3)*(np.cos(0.5*ky*a) + sq3*K*np.sin(0.5*ky*a))) + (np.cos(0.5*ky*a) + sq3*K*np.sin(0.5*ky*a))**2)
    return M2

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

MT1 = (M1(kx,ky)/2)
MT2 = (M2(kx,ky)/2)

# Plot Sigma -1 for K' point
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Matrix Transition with Right Circular Polarized Light', fontsize=15)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.set_xlabel('kx', fontsize=13)
ax.set_ylabel('ky', fontsize=13)
ax.set_zlabel('M^2', fontsize=13)
Ec = ax.plot_surface(kx, ky, MT1, cmap='viridis')
bz = ax.plot(bzx, bzy, 3.5, color='red')
print("Saving matrix transition at K' point as photo")
plt.savefig("hBNMatrixTransitionK'.png")
plt.show()

# Plot Sigma 1 for K point
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Matrix Transition with Left Circular Polarized Light', fontsize=15)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.set_xlabel('kx', fontsize=13)
ax.set_ylabel('ky', fontsize=13)
ax.set_zlabel('M^2', fontsize=13)
Ec = ax.plot_surface(kx, ky, MT2, cmap='viridis')
bz = ax.plot(bzx, bzy, 3.5, color='red')
print("Saving matrix transition at K point as photo")
plt.savefig("hBNMatrixTransitionK.png")
plt.show()