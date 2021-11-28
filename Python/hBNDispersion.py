# -*- coding: utf-8 -*-
"""
Created on Wed May 13 17:56:41 2020

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
t = -1
eps = 0
d = 0.5

# Function
def w(kx,ky):
    W = np.sqrt(1.0+4.0*np.cos(0.5*sq3*kx*a)*np.cos(0.5*ky*a)+4*(np.cos(0.5*ky*a))**2)
    return W
 
def valence(d,kx,ky):
    ev = (-2*t*s*(w(kx,ky)**2) - np.sqrt((d**2) - ((d*s*w(kx,ky))**2) + 4*(t**2)*(w(kx,ky)**2)))/(2*(1 + (s**2)*(w(kx,ky)**2)))
    return ev
 
def conduction(d,kx,ky):
    ec = (-2*t*s*(w(kx,ky)**2) + np.sqrt((d**2) - ((d*s*w(kx,ky))**2) + 4*(t**2)*(w(kx,ky)**2)))/(2*(1 + (s**2)*(w(kx,ky)**2)))
    return ec

# Generate point for 2D point
# K point to Gamma point
kxkg = np.linspace(0, 0, ndata, endpoint=True)
kykg = np.linspace(4*pi/(3*a), 0, ndata, endpoint=True)
Evkg = np.zeros((ndata, 1))
Eckg = np.zeros((ndata, 1))
Evkg = valence(d,kxkg, kykg)
Eckg = conduction(d,kxkg, kykg)

# Gamma point to M point
kxgm = np.linspace(0, 2*pi/sq3, ndata, endpoint=True)
kygm = np.linspace(0, 0, ndata, endpoint=True)
Evgm = np.zeros((ndata, 1))
Ecgm = np.zeros((ndata, 1))
Evgm = valence(d,kxgm, kygm)
Ecgm = conduction(d,kxgm, kygm)

# M point K' point
kxmk = np.linspace(2*pi/sq3, 2*pi/sq3, ndata, endpoint=True)
kymk = np.linspace(0, 2*pi/3, ndata, endpoint=True)
Evmk = np.zeros((ndata, 1))
Ecmk = np.zeros((ndata, 1))
Evmk = valence(d,kxmk, kymk)
Ecmk = conduction(d,kxmk, kymk)

# k as x axis , Evale / Econd / Fermi as y axis
k = np.linspace(0, 2*pi, int(3*ndata), endpoint=True)
Fermi = np.linspace(0, 0, int(3*ndata), endpoint=True)
Evale = np.zeros((int(3*ndata), 1))
Econd = np.zeros((int(3*ndata), 1))
Evale = np.concatenate((Evkg, Evgm, Evmk))
Econd = np.concatenate((Eckg, Ecgm, Ecmk))

# Output 2D Energy Dispersion Data
write = np.column_stack((k.flatten(), Fermi.flatten(), Evale.flatten(), Econd.flatten()))
print('Saving 2D Energy Dispersion data as excel')
np.savetxt('hBNDispersion1D.csv', write, delimiter = ",", header="k,   Fermi Level,   Valence,   Conduction")

# Plot 2D Energy Dispersion
fig, ax = plt.subplots()
plt.title('Energy Dispersion with $\Delta$ = 0.5 eV', fontsize=20)
plt.axis([0, 6.28, -4, 4])
Ev = plt.plot(k, Evale, color='blue')
Ec = plt.plot(k, Econd, color='red')
EF = plt.plot(k, Fermi, color='green', linestyle='--')
ax.set_xticklabels([])
plt.ylabel('Energy (eV)', fontsize=15)
plt.text(0.0, -4.7, '$K$', fontsize=15)
plt.text(2.0, -4.7, '$\Gamma$', fontsize=15)
plt.text(4.1, -4.7, '$M$', fontsize=15)
plt.text(6.2, -4.7, '$K$'"'", fontsize=15)
plt.legend(('Valence Band s = 0', 'Conduction Band s = 0', 'Fermi Level'))
plt.show()
print('Saving 2D Energy Dispersion as photo')
plt.savefig('hBNDispersion1D.png')

# Generate point for 3D point
kmesh = np.linspace(-2*pi, 2*pi, ndata)
kx , ky = np.meshgrid(kmesh, kmesh)

Econd3D = conduction(d,kx,ky)
Evale3D = valence(d,kx,ky)
miu = (d**2 + (d*np.sqrt(d**2 + 4*(t**2)*(w(kx,ky)**2))))/(d**2 + 4*(t**2)*(w(kx,ky)**2) + (d*np.sqrt(d**2 + 4*(t**2)*(w(kx,ky)**2))))

# Output 3D Energy Dispersion Data
write = np.column_stack((kx.flatten(), ky.flatten(), Evale3D.flatten(), Econd3D.flatten(), miu.flatten()))
print('Saving 3D Energy Dispersion data as excel')
np.savetxt('hBNDispersion3D.csv', write, delimiter = ",", header="kx,   ky,   Valence,   Conduction,  miu")

# Plot 3D Energy Dispersion
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Energy Dispersion with $\Delta$ = 0.5 eV', fontsize=20)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel('kx', fontsize=13)
ax.set_ylabel('ky', fontsize=13)
ax.set_zlabel('Energy (eV)', fontsize=13)
Ec = ax.plot_surface(kx, ky, Econd3D)
Ev = ax.plot_surface(kx, ky, Evale3D)
print('Saving 3D Energy Dispersion as photo')
plt.savefig('hBNDispersion3D.png')
plt.show()

# Plot Polarization
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Polarization with $\Delta$ = 0.5 eV', fontsize=20)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel('kx', fontsize=13)
ax.set_ylabel('ky', fontsize=13)
ax.set_zlabel('Energy (eV)', fontsize=13)
P = ax.plot_surface(kx, ky, miu, cmap='coolwarm', vmin=-miu.max(), vmax=miu.max())
print('Saving 3D Energy Dispersion as photo')
plt.savefig('hBNDispersion3D.png')
plt.show()