# -*- coding: utf-8 -*-
"""
Created on Mon May  11 13:22:56 2020

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
s1 = 0.129
t = -2.7
eps = 0

# Function
def w(kx,ky):
    W = np.sqrt(1.0+4.0*np.cos(0.5*sq3*kx*a)*np.cos(0.5*ky*a)+4*(np.cos(0.5*ky*a))**2)
    return W
 
def valence(kx,ky):
    ev = (eps + t*w(kx,ky))/(1.0 + s*w(kx,ky))
    return ev
 
def conduction(kx,ky):
    ec = (eps - t*w(kx,ky))/(1.0 - s*w(kx,ky))
    return ec

def fermi(kx,ky):
    f = 0*kx + 0*ky
    return f

def valence1(kx,ky):
    ev = (eps + t*w(kx,ky))/(1.0 + s1*w(kx,ky))
    return ev
 
def conduction1(kx,ky):
    ec = (eps - t*w(kx,ky))/(1.0 - s1*w(kx,ky))
    return ec

# Generate point for 2D point
# K point to Gamma point
kxkg = np.linspace(0, 0, ndata, endpoint=True)
kykg = np.linspace(4*pi/(3*a), 0, ndata, endpoint=True)
Evkg = np.zeros((ndata, 1))
Eckg = np.zeros((ndata, 1))
Evkg = valence(kxkg, kykg)
Eckg = conduction(kxkg, kykg)
Evkg1 = np.zeros((ndata, 1))
Eckg1 = np.zeros((ndata, 1))
Evkg1 = valence1(kxkg, kykg)
Eckg1 = conduction1(kxkg, kykg)

# Gamma point to M point
kxgm = np.linspace(0, 2*pi/sq3, ndata, endpoint=True)
kygm = np.linspace(0, 0, ndata, endpoint=True)
Evgm = np.zeros((ndata, 1))
Ecgm = np.zeros((ndata, 1))
Evgm = valence(kxgm, kygm)
Ecgm = conduction(kxgm, kygm)
Evgm1 = np.zeros((ndata, 1))
Ecgm1 = np.zeros((ndata, 1))
Evgm1 = valence1(kxgm, kygm)
Ecgm1 = conduction1(kxgm, kygm)

# M point K' point
kxmk = np.linspace(2*pi/sq3, 2*pi/sq3, ndata, endpoint=True)
kymk = np.linspace(0, 2*pi/3, ndata, endpoint=True)
Evmk = np.zeros((ndata, 1))
Ecmk = np.zeros((ndata, 1))
Evmk = valence(kxmk, kymk)
Ecmk = conduction(kxmk, kymk)
Evmk1 = np.zeros((ndata, 1))
Ecmk1 = np.zeros((ndata, 1))
Evmk1 = valence1(kxmk, kymk)
Ecmk1 = conduction1(kxmk, kymk)

# k as x axis , Evale / Econd / Fermi as y axis
k = np.linspace(0, 2*pi, int(3*ndata), endpoint=True)
Fermi = np.linspace(0, 0, int(3*ndata), endpoint=True)
Evale = np.zeros((int(3*ndata), 1))
Econd = np.zeros((int(3*ndata), 1))
Evale = np.concatenate((Evkg, Evgm, Evmk))
Econd = np.concatenate((Eckg, Ecgm, Ecmk))
Evale1 = np.zeros((int(3*ndata), 1))
Econd1 = np.zeros((int(3*ndata), 1))
Evale1 = np.concatenate((Evkg1, Evgm1, Evmk1))
Econd1 = np.concatenate((Eckg1, Ecgm1, Ecmk1))

# Output 2D Energy Dispersion Data
write = np.column_stack((k.flatten(), Fermi.flatten(), Evale.flatten(), Econd.flatten()))
print('Saving 2D Energy Dispersion as excel')
np.savetxt('GrapheneDispersion1D.csv', write, delimiter = ",", header="k,   Fermi Level,   Valence,   Conduction")

# Plot 2D Energy Dispersion
fig, ax = plt.subplots()
plt.title('Energy Dispersion', fontsize=20)
plt.axis([0, 6.28, -9, 15])
Ev = plt.plot(k, Evale, color='blue')
Ec = plt.plot(k, Econd, color='red')
EF = plt.plot(k, Fermi, color='green', linestyle='--')
Ev1 = plt.plot(k, Evale1, color='purple')
Ec1 = plt.plot(k, Econd1, color='orange')
ax.set_xticklabels([])
plt.ylabel('Energy (eV)', fontsize=15)
plt.text(0.0, -11, '$K$', fontsize=15)
plt.text(2.0, -11, '$\Gamma$', fontsize=15)
plt.text(4.1, -11, '$M$', fontsize=15)
plt.text(6.2, -11, '$K$'"'", fontsize=15)
plt.legend(('Valence Band s = 0', 'Conduction Band s = 0', 'Fermi Level', 'Valence Band s = 0.129', 'Conduction Band s = 0.129'))
plt.show()
print('Saving 2D Energy Dispersion as photo')
plt.savefig('GrapheneDispersion1D.png')

# Generate point for 3D point
kmesh = np.linspace(-2*pi, 2*pi, ndata)
kx , ky = np.meshgrid(kmesh, kmesh)

Econd3D = conduction(kx,ky)
Evale3D = valence(kx,ky)
Fermi3D = fermi(kx,ky)

# Output 3D Energy Dispersion Data
write = np.column_stack((kx.flatten(), ky.flatten(), Fermi3D.flatten(), Evale3D.flatten(), Econd3D.flatten()))
print('Saving 3D Energy Dispersion as excel')
np.savetxt('GrapheneDispersion3D.csv', write, delimiter = ",", header="kx,   ky,   Fermi,   Valence,   Conduction")

# Plot 3D Energy Dispersion
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Energy Dispersion', fontsize=20)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel('kx', fontsize=13)
ax.set_ylabel('ky', fontsize=13)
ax.set_zlabel('Energy (eV)', fontsize=13)
Ec = ax.plot_surface(kx, ky, Econd3D, color='red')
Ev = ax.plot_surface(kx, ky, Evale3D, color='blue')
Ef = ax.plot_surface(kx, ky, Fermi3D, color='green')
plt.show()
print('Saving 3D Energy Dispersion as photo')
plt.savefig('GrapheneDispersion3D.png')