# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 18:52:05 2020

@author(s):
    Prof. Herman Batelaan (@unl.edu), Anjaneshwar Ganesan (@unl.edu),
    and Efren A. Serra (@home.working)
"""

import numpy as np
import matplotlib.pyplot as plt

I0 = 1.0
# Mean intensity of light source

elements = 11
# Number of laser outlets

attempts = 10**6
# Number of trials

T = [0,0,1,1,1,1,1,0,0,0,0]
correlated_signal = np.array([ 0. for n in range(elements) ], dtype=np.float)

signal_hist=np.array([0 for n in range(elements)], dtype=np.int)
reference_hist=np.array([0 for n in range(elements)], dtype=np.int)
for k in range(attempts):
    randomized_photon_sources = np.random.randint(0, 11, elements)
    # Randomized photon sources
    randomized_photon_intensities = I0 * np.random.uniform(0., 1., elements)

    signal=np.array([0. for n in range(elements)], dtype=np.float)
    # E_S : signal photon: interacts with the pinhole detector
    reference=np.array([0. for n in range(elements)], dtype=np.float)
    # E_R : reference photon: interacts with the amplitude transmission mask
    for q in randomized_photon_sources:
        signal[q] += randomized_photon_intensities[q]
        reference[q] += randomized_photon_intensities[q]

        signal_hist[q] = signal_hist[q] + 1
        reference_hist[q] = reference_hist[q] + 1

        bucket = np.sum(np.multiply(T, reference))
        # Interaction of photons and transmission mask

        correlated = bucket * np.array(signal)
        # Multiplying the bucket current and the original signal current

        correlated_signal += correlated
        # Sum up correlated signals

plt.figure(1)
plt.plot(correlated_signal/attempts, "ko")
plt.ylabel("Ghost image [Average correlated intensities]")
plt.xlabel("Transmission mask: %s"%T)

plt.figure(2)
plt.bar([n for n in range(elements)], signal_hist/attempts)
plt.ylabel("Signal/pinhole detector photon count")
plt.xlabel("Transmission mask [reversed]: %s"%T[::-1])

plt.figure(3)
plt.bar([n for n in range(elements)], reference_hist/attempts)
plt.ylabel("Reference/bucket detector photon count")
plt.xlabel("Transmission mask: %s"%T)

print("Correlated signal: %s"%correlated_signal)
print("Signal or idler photon intensities: %s"%signal)
print("Reference photon intensities: %s"%reference)