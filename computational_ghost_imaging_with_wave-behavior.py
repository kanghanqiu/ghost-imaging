# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 19:32:52 2020

d = a/N

     |                                 |
     |                                 |
     |                                 +(L,y)
     |                                 |
     |                                 |
     +---------------------------------+(L,d/2)
(0,0)=----------------L----------------+
     +---------------------------------+(L,-d/2)
     |                                 |
     |                                 |
     |                                 |
     |                                 |
@author: Efren A. Serra
"""

from __future__ import print_function

import math
import matplotlib.pyplot as plt
import numpy as np

N : int = 1000
# number of discrete points

D : float = .2 * 10 ** -6
# slit width [micron]
Y : float = 1.0
# Screen width [m]

Dx : float = D / N
Dy : float = Y / N
# deltas

L : int = 100.
# Distance to screen

LAMBDA : float = 650. * 10 ** -9
k : float = (2. * math.pi) / LAMBDA
# Units: nm
# Red lasers typical wavelengths are: 630, 650 and 670 nm
THETA : float = LAMBDA/D

def Intensity(L : float, k : float, x: object, y: object):
    nslits = x.size
    I = np.array([0. for i in range(nslits)])
    for j, yj in np.ndenumerate(y):
        for yk in np.nditer(x):
            delta_kj = math.sqrt(L**2 + (yj-yk)**2)
            I[j] = I[j] + math.cos(k*delta_kj)**2

    return I/nslits

def main():
    x = np.array([i for i in np.arange(-D/2, D/2, Dx)])
    # Mini-slit points
    y = np.array([i for i in np.arange(-1./2, 1./2, Dy)])
    # Screen points

    Iy = Intensity(L, k, x, y)
    plt.title('Single-slit diffraction pattern; No. mini-slits: %d\nSlit widht: 0.2 [micron]'%(N,))
    plt.plot(y[450:550], Iy[450:550], '-x')
    plt.xlabel("y [m]")
    plt.ylabel('Normalized\nIntensity')

if __name__ == "__main__":
    main()