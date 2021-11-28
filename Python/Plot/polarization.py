#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 13:59:13 2021

@author: rudye
"""

# Import module
import numpy as np
import matplotlib.pyplot as plt

# Import datafile .dat / .in to text for pure NbS2
# Python start with index 0, then 1, 2, 3 and so on
xs = np.split(np.loadtxt("./nbs2.degree_pol.dat", unpack=True, usecols=[0]), 1)

sleft1 = np.split(np.loadtxt("./nbs2.degree_pol.dat", unpack=True, usecols=[4]), 1)
sright1 = np.split(np.loadtxt("./nbs2.degree_pol.dat", unpack=True, usecols=[5]), 1)

sleft2 = np.split(np.loadtxt("./nbs2.degree_pol.dat", unpack=True, usecols=[6]), 1)
sright2 = np.split(np.loadtxt("./nbs2.degree_pol.dat", unpack=True, usecols=[7]), 1)

# Interpolation data for pure NbS2
# Changing from text into data type in an array
Xs, Ls1, Rs1, Ls2, Rs2  = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
for i in range(len(xs)):
        Xs = np.append(Xs, xs[i])
        Ls1 = np.append(Ls1, sleft1[i])
        Rs1 = np.append(Rs1, sright1[i])
        Ls2 = np.append(Ls2, sleft2[i])
        Rs2 = np.append(Rs2, sright2[i])
        
xis = np.linspace(Xs.min(), Xs.max(), 30)

# Plotting the optical absorption intensity for pure NbS2
fig, ax = plt.subplots()
plt.axis([0, 0.9603, 0, 1.5])
ax.set_yticks([0, 0.5, 1, 1.5])

plt.plot(xis, Ls1/Ls1.max(), color="red", marker = "v", markersize=12)
plt.plot(xis, Rs1/Rs1.max(), color="red", marker = "^", markersize=12)
plt.plot(xis, Ls2/Ls2.max(), color="blue", marker = "<", markersize=12)
plt.plot(xis, Rs2/Rs2.max(), color="blue", marker = ">", markersize=12)
plt.legend(('LCP $E_L$ = 1.573 eV', 'RCP $E_L$ = 1.573 eV', 'LCP $E_L$ = 2.131 eV', 'RCP $E_L$ = 2.131 eV'), loc=1)
fig.savefig('polarization.eps', format='eps')

