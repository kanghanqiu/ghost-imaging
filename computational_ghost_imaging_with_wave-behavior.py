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
import numexpr as ne

N : int = 100
# number of discrete points

D : float = 1.4e-3
# slit width [micron]
Y : float = 1.
# Screen width [m]

dx : float = D / N
dy : float = Y / N
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
#    E_real = np.array([0. for i in range(nmslits)])
#    E_imag = np.array([1j*0. for i in range(nmslits)])
    wavefront = np.zeros(nmslits, dtype=np.complex128)

    for j, yj in np.ndenumerate(y):
        for yk in np.nditer(x):
            delta_kj = math.sqrt(L ** 2 + (yj-yk) ** 2)
#            E_real[j] += math.cos(k*delta_kj)
#            E_imag[j] += 1j*math.sin(k*delta_kj)
            wavefront[j] += (math.cos(k*delta_kj) + 1j*math.sin(k*delta_kj))

#    return E_real+E_imag
    return wavefront

def Intensity(L : float, k : float, x: object, y: object):
    Ey = E(L, k, x, y)

#    return (Ey * np.conj(Ey))/np.sqrt(np.vdot(Ey, Ey))
#    return ne.evaluate("real(abs(Ey))**2")
    return np.abs(Ey) ** 2

def main():
    x = np.array([i for i in np.arange(-D/2, D/2, dx)])
    # Mini-slit points
    y = np.array([i for i in np.arange(-Y/2, Y/2, dy)])
    # Screen points

    Iy = Intensity(L, k, x, y)
    plt.title('Single-slit diffraction pattern; no. mini-slits: %d\nSlit width: %.2e [m]; $\lambda$: %e [m]'%(N,D,LAMBDA))
    plt.plot(y, Iy, '-x')
    plt.xlabel("y [m]")
    plt.ylabel('Normalized\nIntensity')

if __name__ == "__main__":
    main()