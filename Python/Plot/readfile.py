#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 13:59:13 2021

@author: rudye
"""

# Import module
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Import datafile .dat / .in for reading to text for pure NbS2
# Python start with index 0, then 1, 2, 3 and so on
xs = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[0]), 1)
ys = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[1]), 1)

sleft1 = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[3]), 1)
sright1 = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[4]), 1)

sleft2 = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[5]), 1)
sright2 = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[6]), 1)

# Interpolation data for pure NbS2
# Changing from text into data type in an array
Xs, Ys, Ls1, Rs1, Ls2, Rs2  = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
for i in range(len(xs)):
        Xs = np.append(Xs, xs[i])
        Ys = np.append(Ys, ys[i])
        Ls1 = np.append(Ls1, sleft1[i])
        Rs1 = np.append(Rs1, sright1[i])
        Ls2 = np.append(Ls2, sleft2[i])
        Rs2 = np.append(Rs2, sright2[i])

# Interpolate to 1000 data
xis = np.linspace(Xs.min(), Xs.max(), 1000)
yis = np.linspace(Ys.min(), Ys.max(), 1000)

# For making a grid mesh for 3D plot
zis1 = griddata((Ys, Xs), Ls1/Ls1.max(), (xis[None,:], yis [:,None]), method='cubic')
zis2 = griddata((Ys, Xs), Rs1/Rs1.max(), (xis[None,:], yis [:,None]), method='cubic')
zis3 = griddata((Ys, Xs), Ls2/Ls2.max(), (xis[None,:], yis [:,None]), method='cubic')
zis4 = griddata((Ys, Xs), Rs2/Rs2.max(), (xis[None,:], yis [:,None]), method='cubic')

# Plotting the matrix element transition for pure NbS2 with heatmap
fig, sx2 = plt.subplots()
sx2.axis([-0.9, 0.9, -0.9, 0.9])
sx2.set_xticklabels([])
sx2.set_yticklabels([])
sx2.text(-0.02, -0.02, '$\Gamma$', color="white", fontsize=15)
sx2.text(-0.67, 0.45, "K'", color="white", fontsize=15)
sx2.text(0.67, -0.45, "K", color="white", fontsize=15)
sx2.set_xticklabels([])
sx2.set_yticklabels([])
s_nondop_right1 = sx2.contourf(xis, yis, zis2, 80, cmap=plt.cm.gnuplot, vmin=0, vmax=1)
fig.colorbar(s_nondop_right1, ticks = np.linspace(0, 1, 5), ax=sx2)
fig.savefig('NbS2_pure_right1.eps', format='eps')

fig, sx1 = plt.subplots()
sx1.axis([-0.9, 0.9, -0.9, 0.9])
sx1.set_xticklabels([])
sx1.set_yticklabels([])
sx1.text(-0.02, -0.02, '$\Gamma$', color="white", fontsize=15)
sx1.text(-0.67, 0.45, "K'", color="white", fontsize=15)
sx1.text(0.67, -0.45, "K", color="white", fontsize=15)
sx1.set_xticklabels([])
sx1.set_yticklabels([])
s_nondop_left1 = sx1.contourf(xis, yis, zis1, 80, cmap=plt.cm.gnuplot, vmin=0, vmax=1)
fig.colorbar(s_nondop_right1, ticks = np.linspace(0, 1, 5), ax=sx1)
fig.savefig('NbS2_pure_left1.eps', format='eps')

fig, sx3 = plt.subplots()
sx3.axis([-0.9, 0.9, -0.9, 0.9])
sx3.set_xticklabels([])
sx3.set_yticklabels([])
sx3.text(-0.02, -0.02, '$\Gamma$', color="white", fontsize=15)
sx3.text(-0.67, 0.45, "K'", color="white", fontsize=15)
sx3.text(0.67, -0.45, "K", color="white", fontsize=15)
sx3.set_xticklabels([])
sx3.set_yticklabels([])
s_nondop_left2 = sx3.contourf(xis, yis, zis3, 80, cmap=plt.cm.gnuplot, vmin=0, vmax=1)
fig.colorbar(s_nondop_right1, ticks = np.linspace(0, 1, 5), ax=sx3)
fig.savefig('NbS2_pure_left2.eps', format='eps')

fig, sx4 = plt.subplots()
sx4.axis([-0.9, 0.9, -0.9, 0.9])
sx4.set_xticklabels([])
sx4.set_yticklabels([])
sx4.text(-0.02, -0.02, '$\Gamma$', color="white", fontsize=15)
sx4.text(-0.67, 0.45, "K'", color="white", fontsize=15)
sx4.text(0.67, -0.45, "K", color="white", fontsize=15)
sx4.set_xticklabels([])
sx4.set_yticklabels([])
s_nondop_right2 = sx4.contourf(xis, yis, zis4, 80, cmap=plt.cm.gnuplot, vmin=0, vmax=1)
fig.colorbar(s_nondop_right1, ticks = np.linspace(0, 1, 5), ax=sx4)
fig.savefig('NbS2_pure_right2.eps', format='eps')
