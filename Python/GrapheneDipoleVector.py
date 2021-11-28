# -*- coding: utf-8 -*-
"""
Created on Sun May 10 12:38:52 2020

@author: rudye layton
"""

import numpy as np
import matplotlib.pyplot as mpl

# Input 
ndata = 200
pi = np.pi
sq3 = np.sqrt(3)
a = 1
s = 0
t = -2.7
eps = 0

# Function
def w(kx,ky):
    W = np.sqrt(1.0 + 4.0*np.cos(0.5*sq3*kx*a)*np.cos(0.5*ky*a) + 4*(np.cos(0.5*ky*a))**2)
    return W

def Dx(kx,ky):
    Dx = (np.cos(0.5*ky*a)*np.cos(sq3*0.5*kx*a) - np.cos(ky*a))
    return Dx

def Dy(kx,ky):
    Dy = (sq3*np.sin(0.5*ky*a)*np.sin(sq3*0.5*kx*a))
    return Dy

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

# Generate point
kmesh = np.linspace(-2*pi, 2*pi, ndata)
kx , ky = np.meshgrid(kmesh, kmesh)
k = np.sqrt(kx**2 + ky**2)

# Vector Field
Vx, Vy = -Dx(kx,ky)/w(kx,ky), -Dy(kx,ky)/w(kx,ky)

# Output Data
write = np.column_stack((kx.flatten(), ky.flatten(), Vx.flatten(), Vy.flatten()))
print('Saving dipole vector as excel')
np.savetxt('GrapheneDipoleVector.csv', write, delimiter = ",", header="kx,   ky,   Vx,   Vy")

# Plot
fig, ax = mpl.subplots()
mpl.axis([-7, 7, -7, 7])
mpl.title('Dipole Vector', fontsize=20)
ax.set_xticklabels([])
ax.set_yticklabels([])
mpl.text(-0.2, -0.2, '$\Gamma$', fontsize=15)
mpl.text(-0.32, -8, 'kx', fontsize=15)
mpl.text(-8, -0.25, 'ky', fontsize=15)
bz = mpl.plot(bzx, bzy, color='red')
vf = mpl.streamplot(kx, ky, Vx, Vy, arrowsize=1.5)
mpl.show()
print('Saving dipole vector as photo')
mpl.savefig('GrapheneDipoleVector.png')