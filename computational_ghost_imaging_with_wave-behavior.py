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

D : float = 5. * 10 ** -6
# slit width
N : int = 100
# number of discrete points inside slit
Dx : float = D / N
Dy : float = 1./ N
# Units: micron
L : int = 100.
# Units: meter
LAMBDA : float = 650. * 10 ** -9
k : float = (2. * math.pi) / LAMBDA
# Units: nm
# Red lasers typical wavelengths are: 630, 650 and 670 nm
THETA : float = LAMBDA/D

def Intensity(L : float, k : float, x: object, y: object):
    I = np.array([0. for i in range(x.size)])
    for j, yj in np.ndenumerate(y):
        for yk in np.nditer(x):
            delta_kj = math.sqrt(L**2 + (yj-yk)**2)
            I[j] = I[j] + math.cos(k*delta_kj)**2

    return I

def main():
    x = np.array([i for i in np.arange(-D/2, D/2, Dx)])
    # Mini-slit points
    y = np.array([i for i in np.arange(-1./2, 1./2, Dy)])
    # Screen points

    Iy = Intensity(L, k, x, y)
    plt.title('Single-slit diffraction pattern')
    plt.plot(y, Iy, '-x')
    plt.xlabel("y [m]")
    plt.ylabel('Normalized\nIntensity')

if __name__ == "__main__":
    main()