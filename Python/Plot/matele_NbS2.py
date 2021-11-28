#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 13:59:13 2021

@author: rudye
"""

# This program is plotting matrix element transition for pure NbX2
# Uncomment part that need to use for NbSe2 or NbTe2
# Brillouin Zone is still manual , so please change all the exact coordinate based on the material

# Import module
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Import datafile .dat / .in to text for pure NbS2
xs = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[0]), 1)
ys = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[1]), 1)

sleft1 = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[3]), 1)
sright1 = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[4]), 1)

sleft2 = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[5]), 1)
sright2 = np.split(np.loadtxt("./NbS2.mopt_ALL.dat", unpack=True, usecols=[6]), 1)

## Import datafile .dat / .in to text for doped NbS2
#xsd = np.split(np.loadtxt("./NbS2.doped_mopt_ALL.dat", unpack=True, usecols=[0]), 1)
#ysd = np.split(np.loadtxt("./NbS2.doped_mopt_ALL.dat", unpack=True, usecols=[1]), 1)
#
#sdleft1 = np.split(np.loadtxt("./NbS2.doped_mopt_ALL.dat", unpack=True, usecols=[3]), 1)
#sdright1 = np.split(np.loadtxt("./NbS2.doped_mopt_ALL.dat", unpack=True, usecols=[4]), 1)

#Interpolation for pure NbS2
Xs, Ys, Ls1, Rs1, Ls2, Rs2  = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
for i in range(len(xs)):
        Xs = np.append(Xs, xs[i])
        Ys = np.append(Ys, ys[i])
        Ls1 = np.append(Ls1, sleft1[i])
        Rs1 = np.append(Rs1, sright1[i])
        Ls2 = np.append(Ls2, sleft2[i])
        Rs2 = np.append(Rs2, sright2[i])
        
xis = np.linspace(Xs.min(), Xs.max(), 1000)
yis = np.linspace(Ys.min(), Ys.max(), 1000)

zis1 = griddata((Ys, Xs), Ls1/Ls1.max(), (xis[None,:], yis [:,None]), method='cubic')
zis2 = griddata((Ys, Xs), Rs1/Rs1.max(), (xis[None,:], yis [:,None]), method='cubic')
zis3 = griddata((Ys, Xs), Ls2/Ls2.max(), (xis[None,:], yis [:,None]), method='cubic')
zis4 = griddata((Ys, Xs), Rs2/Rs2.max(), (xis[None,:], yis [:,None]), method='cubic')

##Interpolation for doped NbS2
#Xsd, Ysd, Lsd1, Rsd1 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
#for i in range(len(xs)):
#        Xsd = np.append(Xsd, xs[i])
#        Ysd = np.append(Ysd, ys[i])
#        Lsd1 = np.append(Lsd1, sdleft1[i])
#        Rsd1 = np.append(Rs1, sdright1[i])
#        
#xisd = np.linspace(Xsd.min(), Xsd.max(), 1000)
#yisd = np.linspace(Ysd.min(), Ysd.max(), 1000)
#
#zisd1 = griddata((Ysd, Xsd), Lsd1/Lsd1.max(), (xisd[None,:], yisd [:,None]), method='cubic')
#zisd2 = griddata((Ysd, Xsd), Rsd1/Rsd1.max(), (xisd[None,:], yisd [:,None]), method='cubic')

# Brillouin Zone pure NbS2
kbzxs1 = np.linspace(-0.5677, -0.5677, 1000)
kbzys1 = np.linspace(0.3311, -0.298, 1000)

kbzxs2 = np.linspace(0.5735, 0.5735, 1000)
kbzys2 = np.linspace(-0.298, 0.3311, 1000)

kbzxs3 = np.linspace(0.5735, 0.0131, 1000)
kbzys3 = np.linspace(0.3311, 0.6623, 1000)

kbzxs4 = np.linspace(0.0131, -0.5677, 1000)
kbzys4 = np.linspace(0.6623, 0.3311, 1000)

kbzxs5 = np.linspace(-0.5677, 0.0131, 1000)
kbzys5 = np.linspace(-0.298, -0.6292, 1000)

kbzxs6 = np.linspace(0.0131, 0.5735, 1000)
kbzys6 = np.linspace(-0.6292, -0.298, 1000)

bzxs = np.concatenate((kbzxs1,kbzxs5,kbzxs6,kbzxs2,kbzxs3,kbzxs4))
bzys = np.concatenate((kbzys1,kbzys5,kbzys6,kbzys2,kbzys3,kbzys4))

## Brillouin Zone doped NbS2
#kbzxsd1 = np.linspace(-0.5677, -0.5677, 1000)
#kbzysd1 = np.linspace(0.3311, -0.298, 1000)
#
#kbzxsd2 = np.linspace(0.5735, 0.5735, 1000)
#kbzysd2 = np.linspace(-0.298, 0.3311, 1000)
#
#kbzxsd3 = np.linspace(0.5735, 0.0131, 1000)
#kbzysd3 = np.linspace(0.3311, 0.6623, 1000)
#
#kbzxsd4 = np.linspace(0.0131, -0.5677, 1000)
#kbzysd4 = np.linspace(0.6623, 0.3311, 1000)
#
#kbzxsd5 = np.linspace(-0.5677, 0.0131, 1000)
#kbzysd5 = np.linspace(-0.298, -0.6292, 1000)
#
#kbzxsd6 = np.linspace(0.0131, 0.5735, 1000)
#kbzysd6 = np.linspace(-0.6292, -0.298, 1000)
#
#bzxsd = np.concatenate((kbzxsd1,kbzxsd5,kbzxsd6,kbzxsd2,kbzxsd3,kbzxsd4))
#bzysd = np.concatenate((kbzysd1,kbzysd5,kbzysd6,kbzysd2,kbzysd3,kbzysd4))

# Plotting the matrix element transition for pure NbS2
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
bz2 = sx2.plot(bzxs, bzys, linewidth=0.2, color='white')
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
bz1 = sx1.plot(bzxs, bzys, linewidth=0.2, color='white')
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
bz3 = sx3.plot(bzxs, bzys, linewidth=0.2, color='white')
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
bz4 = sx4.plot(bzxs, bzys, linewidth=0.2, color='white')
fig.colorbar(s_nondop_right1, ticks = np.linspace(0, 1, 5), ax=sx4)
fig.savefig('NbS2_pure_right2.eps', format='eps')

## Plotting the matrix element transition for doped NbS2
#fig, sdx2 = plt.subplots()
#sdx2.axis([-0.9, 0.9, -0.9, 0.9])
#sdx2.set_xticklabels([])
#sdx2.set_yticklabels([])
#sdx2.text(-0.05, -0.05, '$\Gamma$', color="white", fontsize=15)
#sdx2.text(-0.67, 0.45, "K'", color="white", fontsize=15)
#sdx2.text(0.67, -0.45, "K", color="white", fontsize=15)
#sdx2.set_xticklabels([])
#sdx2.set_yticklabels([])
#sd_doped_right1 = sdx2.contourf(xisd, yisd, zisd2, 100, cmap=plt.cm.gnuplot, vmin=0, vmax=1)
#bzd2 = sdx2.plot(bzxsd, bzysd, linewidth=0.2, color='white')
#fig.colorbar(sd_doped_right1, ticks = np.linspace(0, 1, 5), ax=sdx2)
#fig.savefig('NbS2_doped_right1.eps', format='eps')
#
#fig, sdx1 = plt.subplots()
#sdx1.axis([-0.9, 0.9, -0.9, 0.9])
#sdx1.set_xticklabels([])
#sdx1.set_yticklabels([])
#sdx1.text(-0.05, -0.05, '$\Gamma$', color="white", fontsize=15)
#sdx1.text(-0.67, 0.45, "K'", color="white", fontsize=15)
#sdx1.text(0.67, -0.45, "K", color="white", fontsize=15)
#sdx1.set_xticklabels([])
#sdx1.set_yticklabels([])
#sd_doped_left1 = sdx1.contourf(xisd, yisd, zisd1, 100, cmap=plt.cm.gnuplot, vmin=0, vmax=1)
#bzd1 = sdx1.plot(bzxsd, bzysd, linewidth=0.2, color='white')
#fig.colorbar(sd_doped_left1, ticks = np.linspace(0, 1, 5), ax=sdx1)
#fig.savefig('NbS2_doped_left1.eps', format='eps')