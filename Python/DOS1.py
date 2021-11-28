# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

t = 1
eta = 0.02                      # eta kecil nk nw besar
eps0 = 0
krange = [-np.pi,np.pi]
wrange = [-6,6]
nk = 600                        # jumlah titik
nw = 300

kk = np.linspace(krange[0] , krange[1], num=nk)
ww = np.linspace(wrange[0] , wrange[1], num=nw)

Green = np.zeros(shape=(len(ww) , len(kk)) , dtype=np.complex64)
for iw in range(len(ww)):
     for ik in range(len(kk)):
         H = -2 * t *np.cos(kk[ik])
         Green[iw,ik] = 1 / (complex(ww[iw] - H,eta))
         
sumk = np.sum(Green,axis=1)
DOS = - sumk.imag / np.pi / nk

plt.plot(ww,DOS)
plt.xlabel('$\omega$ (eV)')
plt.ylabel("DOS $\omega$")