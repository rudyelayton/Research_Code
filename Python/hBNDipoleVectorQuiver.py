# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 16:01:11 2020

@author: ASUS
"""

import numpy as np
import matplotlib.pyplot as mpl

# Input 
ndata = 200
pi = np.pi
sq3 = np.sqrt(3)
a = 1
s = 0
t = -1.1
eps = 0
d = 2
j = 1j

#R Position
R1x = a/sq3
R1y = 0
R2x = -0.5*a/sq3
R2y = 0.5*a
R3x = -0.5*a/sq3
R3y = -0.5*a

# Function
def w(kx,ky):
    W = np.sqrt(1.0 + 4.0*np.cos(0.5*sq3*kx*a)*np.cos(0.5*ky*a) + 4*(np.cos(0.5*ky*a))**2)
    return W

def fconjugate(kx,ky):
    fconjugate = (np.exp(-j*kx*a/sq3)) + (np.exp(j*kx*a*0.5/sq3))*(2*np.cos(ky*a*0.5))
    return fconjugate

def X(kx,ky):
    x = d/2 + np.sqrt(((d**2)/4) + (t**2)*(w(kx,ky)**2))
    return x

def sume1x(kx,ky):
    sume1x = (np.exp(j*kx*R1x)*np.exp(j*ky*R1y))*R1x + (np.exp(j*kx*R2x)*np.exp(j*ky*R2y))*R2x + (np.exp(j*kx*R3x)*np.exp(j*ky*R3y))*R3x 
    return sume1x

def sume1y(kx,ky):
    sume1y = (np.exp(j*kx*R1x)*np.exp(j*ky*R1y))*R1y + (np.exp(j*kx*R2x)*np.exp(j*ky*R2y))*R2y + (np.exp(j*kx*R3x)*np.exp(j*ky*R3y))*R3y
    return sume1y

def sume2x(kx,ky):
    sume2x = (np.exp(-j*kx*R1x)*np.exp(-j*ky*R1y))*R1x + (np.exp(-j*kx*R2x)*np.exp(-j*ky*R2y))*R2x + (np.exp(-j*kx*R3x)*np.exp(-j*ky*R3y))*R3x 
    return sume2x

def sume2y(kx,ky):
    sume2y = (np.exp(-j*kx*R1x)*np.exp(-j*ky*R1y))*R1y + (np.exp(-j*kx*R2x)*np.exp(-j*ky*R2y))*R2y + (np.exp(-j*kx*R3x)*np.exp(-j*ky*R3y))*R3y
    return sume2y

def Dx(kx,ky):
    Dx = (t**2)*((fconjugate(kx,ky))**2)*sume1x(kx,ky) + ((X(kx,ky))**2)*sume2x(kx,ky) 
    return Dx

def Dy(kx,ky):
    Dy = (t**2)*((fconjugate(kx,ky))**2)*sume1y(kx,ky) + ((X(kx,ky))**2)*sume2y(kx,ky)
    return Dy

def Dxf(kx,ky):
    Dxf = (-sq3/a)*(1/((t**2)*((w(kx,ky))**2) + (X(kx,ky)**2)))*(Dx(kx,ky))
    return Dxf

def Dyf(kx,ky):
    Dyf = (-sq3/a)*(1/((t**2)*((w(kx,ky))**2) + (X(kx,ky)**2)))*(Dy(kx,ky))
    return Dyf

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
Vx, Vy = Dxf(kx,ky), Dyf(kx,ky)

# Plot
fig, ax = mpl.subplots()
mpl.axis([-7, 7, -7, 7])
mpl.title('hBN Dipole Vector with $\Delta$ = 2', fontsize=20)
ax.set_xticklabels([])
ax.set_yticklabels([])
mpl.text(-0.2, -0.2, '$\Gamma$', fontsize=15)
mpl.text(-0.32, -8, 'kx', fontsize=15)
mpl.text(-8, -0.25, 'ky', fontsize=15)
bz = mpl.plot(bzx, bzy, color='black')
vfr = mpl.quiver(kx, ky, Vx.real, Vy.real, color='blue')
vfi = mpl.quiver(kx, ky, Vx.imag, Vy.imag, color='green')
mpl.show()