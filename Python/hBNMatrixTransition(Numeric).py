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
t = -1.1
eps = 0
d = 0.5
j = 1j 
sigma = 1

#R Position
R1x = a/sq3
R1y = 0
R2x = -0.5*a/sq3
R2y = 0.5*a
R3x = -0.5*a/sq3
R3y = -0.5*a

# Function
def w(kx,ky):
    w = np.sqrt(1.0 + 4.0*np.cos(0.5*sq3*kx*a)*np.cos(0.5*ky*a) + 4*(np.cos(0.5*ky*a))**2)
    return w

def fconjugate(kx,ky):
    fconjugate = (np.exp(-j*kx*a/sq3)) + (np.exp(j*kx*a*0.5/sq3))*(2*np.cos(ky*a*0.5))
    return fconjugate

def X(kx,ky):
    x = d/2 + np.sqrt(((d**2)/4) + (t**2)*(w(kx,ky)**2))
    return x

def sume1(kx,ky):
    sume1 = np.dot(1, ((np.exp(j*kx*R1x)*np.exp(j*ky*R1y))*R1x + (np.exp(j*kx*R2x)*np.exp(j*ky*R2y))*R2x + (np.exp(j*kx*R3x)*np.exp(j*ky*R3y))*R3x)) + np.dot(-j*sigma, ((np.exp(j*kx*R1x)*np.exp(j*ky*R1y))*R1y + (np.exp(j*kx*R2x)*np.exp(j*ky*R2y))*R2y + (np.exp(j*kx*R3x)*np.exp(j*ky*R3y))*R3y))
    return sume1

def sume2(kx,ky):
    sume2 = np.dot(1, ((np.exp(-j*kx*R1x)*np.exp(-j*ky*R1y))*R1x + (np.exp(-j*kx*R2x)*np.exp(-j*ky*R2y))*R2x + (np.exp(-j*kx*R3x)*np.exp(-j*ky*R3y))*R3x)) + np.dot(-j*sigma, ((np.exp(-j*kx*R1x)*np.exp(-j*ky*R1y))*R1y + (np.exp(-j*kx*R2x)*np.exp(-j*ky*R2y))*R2y + (np.exp(-j*kx*R3x)*np.exp(-j*ky*R3y))*R3y))
    return sume2

def M(kx,ky):
    M = (t**2)*((fconjugate(kx,ky))**2)*sume1(kx,ky) + (X(kx,ky)**2)*sume2(kx,ky)
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


Mf = (-sq3/a)*(1/((t**2)*((w(kx,ky))**2) + (X(kx,ky)**2)))*M(kx,ky)
Mfconj = np.conj(Mf)

MT = np.abs(Mf*Mfconj)

# Plot Sigma 1 for K point
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Matrix Transition LCP $\Delta$ = 0.5', fontsize=15)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.set_xlabel('kx', fontsize=13)
ax.set_ylabel('ky', fontsize=13)
ax.set_zlabel('M^2', fontsize=13)
Ec = ax.plot_surface(kx, ky, MT.real, cmap='viridis')
bz = ax.plot(bzx, bzy, 10, color='red')
cb = fig.colorbar(Ec)
print("Saving matrix transition at K point as photo")
plt.savefig("hBNMatrixTransitionK(Numeric).png")
plt.show()