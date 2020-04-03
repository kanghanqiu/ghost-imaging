# -*- coding: utf-8 -*-
"""

                                     +
                                     +
                                     +  pinhole detector
                                     +
         \                           +     ^    (scanning)
  E(rho,t)\  Es(rho,t)     E1(rho,t) + D---|-----+
   ======> \  ======>        ======> +     |     |
            \<---------------------->+     |     |
           ^                         +     v     |
 Er(rho,t) |                         +           |
           |                                     |
           |                                     v
           |                                +----------+
 E2(rho,t) |    object, T(rho)              |correlator|
           v       /                        +----------+
      +=========+ <+                              ^
      +---------+                                 |
      |         |                                 |
      \         /                                 |
       +-------+                                  |
           +                                      |
           |                                      |
           +--------------------------------------+
Created on Fri Mar 20 19:32:52 2020

@author: Efren A. Serra
"""

from __future__ import print_function

import math
import matplotlib.pyplot as plt
import numpy as np
import numexpr as ne

N : int = 2 ** 8
# number of discrete points

D : float = 20e-6
# slit width [micron]

Y : float = 1.0e-2
# Screen width [m]
W : float = 2 * D

drho_s : float = D / N
drho_1 : float = Y / N
drho_2 : float = W / N
# Deltas

z_direct : int = 1.0
# Distance to pinhole detector, L [m]

wavelength : float = 650.e-9
# Units: nm
# Red lasers typical wavelengths are: 630, 650 and 670 nm
dtheta : float = wavelength/D

def K(L : float, k : float, x: object, y: object):
    N = x.size
#    E_real = np.array([0. for i in range(nmslits)])
#    E_imag = np.array([1j*0. for i in range(nmslits)])
    wavefront = np.zeros(N, dtype=np.complex128)

    for j, yj in np.ndenumerate(y):
        for yk in np.nditer(x):
            delta_kj = math.sqrt(L ** 2 + (yj-yk) ** 2)
#            E_real[j] += math.cos(k*delta_kj)
#            E_imag[j] += 1j*math.sin(k*delta_kj)
            wavefront[j] += (math.cos(k*delta_kj) + 1.0j*math.sin(k*delta_kj))

#    return E_real+E_imag
    return wavefront

def Intensity(L : float, k : float, x: object, y: object):
    Ey = K(L, k, x, y)

#    return (Ey * np.conj(Ey))/np.sqrt(np.vdot(Ey, Ey))
#    return ne.evaluate("real(abs(Ey))**2")
    return np.abs(Ey) ** 2 / np.sum(np.abs(Ey) ** 2)

def E(N, wavelength):
    """
    This represents a complex-valued uniform random process.

    Parameters
    ----------
        N : int
        Number of random sources for cw laser, e.g., Es(rho, t) field.
        
    Returns
    -------
    """
    k0 = 2 * np.pi / wavelength
    return np.random.uniform(0., 1., N)*np.exp(1.0j * k0)

def propagate_direct(Em, rho, rho_p, L, wavelength):
    """
    
    Parameters
    ----------
    Em (m = S or R)
    rho
    rho_p
    L
    wavelength:
    Returns
    -------
    """
    k0 = 2 * np.pi / wavelength
    El = np.zeros(N, dtype=np.complex128)

    for j, rhoj in np.ndenumerate(rho):
        for k, rho_pk in np.ndenumerate(rho_p):
            quadphase_1st = np.exp(1.0j * k0 * (L + (rhoj-rho_pk) ** 2 / (2 * L))) / (1.0j * 2 * np.pi * L)
            # eq. (1)
            El[j] += (Em[k] * quadphase_1st)

    return El

def main():
    n_trials = 10 ** 1
    rho_s = np.array([i for i in np.arange(-D/2, D/2, drho_s)])
    rho_1 = np.array([i for i in np.arange(-Y/2, Y/2, drho_1)])
    rho_2 = np.array([i for i in np.arange(-W/2, W/2, drho_2)])

    T_rho = np.where(rho_2 <= -(W / 4.0), 0.0, 1.0)
    T_rho = np.where(rho_2 >= (W / 4.0), 0.0, 1.0 * T_rho)

    I2 = np.zeros(N)
    C_rho1 = np.zeros(N, dtype=np.complex128)
    for j in range(N):
        for n in range(n_trials):
            Es = E(N, wavelength)
            E1 = propagate_direct(Es / math.sqrt(2), rho_s, rho_1, z_direct, wavelength)
            I1 = E1 * np.conj(E1)
#            I1 = ne.evaluate("real(abs(E1)**2)")

            E2 = propagate_direct(Es / math.sqrt(2), rho_s, rho_2, z_direct, wavelength)
            I2 = E2 * np.conj(E2)
#            I2 = ne.evaluate("real(abs(E2)**2)")
    
            # Multiplying the bucket current and the original signal current
            C_rho1[j] += (np.sum(I2 * T_rho) * I1[j])

    plt.figure(1)
    plt.title('Intensity pattern, bucket detector, center $\\rho_{1}$.\n$\lambda$: %e [m]'%(wavelength))
    plt.plot(rho_1, C_rho1 / np.sqrt(np.vdot(C_rho1, C_rho1)) / n_trials, '-')
    plt.xlabel("$\\rho_{1}$ [m]")
    plt.ylabel('Intensity [normalized units]')

if __name__ == "__main__":
    main()