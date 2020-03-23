# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 19:32:52 2020


      |                                 |
      |                                 |
      |                                 +(L,y)
      |                                 |
      |                                 |
   +  +---------------------------------+(L,d/2)
 D |  +----------------L----------------+(L,0)
   +  +---------------------------------+(L,-d/2)
      |                                 |
      |                                 |
      |                                 |
      |                                 |

d = D/N

@author: Efren A. Serra
"""

from __future__ import print_function

import math
import matplotlib.pyplot as plt
import numpy as np

N : int = 100
# number of discrete points

D : float = 6.e-4
# slit width [micron]
Y : float = 1.
# Screen width [m]

Dx : float = D / N
Dy : float = Y / N
# deltas

L : int = 100.
# Distance to screen

LAMBDA : float = 650.e-9
k : float = (2. * math.pi) / LAMBDA
# Units: nm
# Red lasers typical wavelengths are: 630, 650 and 670 nm
THETA : float = LAMBDA/D

def E(L : float, k : float, x: object, y: object):
    nmslits = x.size
    E_real = np.array([0. for i in range(nmslits)])
    E_imag = np.array([1j*0. for i in range(nmslits)])

    for j, yj in np.ndenumerate(y):
        for yk in np.nditer(x):
            delta_kj = math.sqrt(L ** 2 + (yj-yk) ** 2)
            E_real[j] += math.cos(k*delta_kj)
            E_imag[j] += 1j*math.sin(k*delta_kj)

    return E_real+E_imag

def Intensity(L : float, k : float, x: object, y: object):
    Ey = E(L, k, x, y)

    return (Ey * np.conj(Ey))/np.sqrt(np.vdot(Ey, Ey))

def main():
    x = np.array([i for i in np.arange(-D/2, D/2, Dx)])
    # Mini-slit points
    y = np.array([i for i in np.arange(-Y/2, Y/2, Dy)])
    # Screen points

    Iy = Intensity(L, k, x, y)
    plt.title('Single-slit diffraction pattern; no. mini-slits: %d\nSlit width: %.2e [m]; $\lambda$: %e [m]'%(N,D,LAMBDA))
    plt.plot(y, Iy, '-x')
    plt.xlabel("y [m]")
    plt.ylabel('Normalized\nIntensity')

if __name__ == "__main__":
    main()