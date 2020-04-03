# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 14:02:41 2020

@author: Efren A. Serra
"""
import numpy as np
import numpy.fft as np_fft
import matplotlib.pyplot as plt
#class Wavefront():
I=1000
#% 1000 Watts per square meter = 1 Watt/mm (per meter of slit length).

N=2 ** 12
#% Number of sample points.

L=0.66e-6
#% Wavelength (red light).

D=1000e-6
#% Total x distance sampled (1mm).

W=20e-6
#% Slit width (20um).

dx=D/N
#% The x-domain sample interval (approx).

dq=1/D
#% The q-domain sample interval (q = theta/lambda).

dtheta=L*dq
#% Angular displacement interval.

x=np.zeros(N)
#% Construct x vector proportional to the field amplitude.

W2=np.round(W/2/dx)
x[int(N/2-W2) : int(N/2+W2)]=np.sqrt(I)
#% Set field strength constant over the slit width.

xf=dx*np_fft.fftshift(np_fft.fft(x))
#% FFT scaled and shifted - Transform domain is q = theta/lambda).

x_domain_total_power=dx*np.sum(np.abs(x) ** 2)
#% Compare the total power in x-domain with that

q_domain_total_power=dq*np.sum(np.abs(xf) ** 2)
#% of the q-domain to confirm scaling is correct.
plt.subplot(2,1,1)
plt.plot(np.array([i for i in np.arange(int(-N/2),int(N/2))])*dx*1000, np.abs(x) ** 2/1000)
#% Plot the x-domain intensity (converting to per mm).
plt.xlabel('Slit Position (mm)')
plt.ylabel('Watts/mm (per m length)')
plt.title ('Slit power density versus position : Width=20um : Lambda=0.66um.')

plt.subplot(212)
plt.plot(np.array([i for i in np.arange(int(-N/2),int(N/2))])*dtheta,np.abs(xf) ** 2/L)
#% Plot q-domain intensity (converting to per radian).
plt.xlabel('Angular Position (radians)')
plt.ylabel('Watts/rad (per m length)')
plt.title ('Diffraction pattern, power density versus angular position.')

"""
Source https://www.physicsforums.com/threads/fft-single-slit-diffraction.377324/
"""