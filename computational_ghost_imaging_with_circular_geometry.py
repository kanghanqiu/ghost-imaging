# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 18:52:05 2020

@author(s):
    Prof. Herman Batelaan (@unl.edu), Anjaneshwar Ganesan (@unl.edu),
    and Efren A. Serra (@home.work)
"""

import math
import numpy as np
import matplotlib.pyplot as plt

I0 = 1.0
# Mean intensity of light source

attempts = 10**5
# Number of trials

def object_T(theta: float, a: float, b: float) -> float:
    """
    The object transmission mask lying in interval: [a, b]
    on the unit circle.
    """
    if a <= theta and theta <= b:
        return 1
    else:
        return 0

def pinhole_number(theta: float, pinhole_spacing: float=0.5) -> int:
    """
    The pinhole detector position between 0 and 180 degrees, every .5 degrees,
    on the unit circle.
    """
    return math.floor(theta / pinhole_spacing)

n_pinholes = int(180./.5)
# One pinhole detector every 1/2 degrees
object_T_mask_range = [ (225., 255.), (255., 285.), (285., 315.) ]

bucket_n = [(n_pinholes-(n+1)) for n in range(n_pinholes)]

for i, T_mask_range in enumerate(object_T_mask_range):
    T = np.array([object_T(theta, *T_mask_range) for theta in np.arange(180.0, 360., .5)])
    correlated_signal = np.array([ 0. for n in range(n_pinholes) ], dtype=np.float)
    signal_hist=np.array([0 for n in range(n_pinholes)], dtype=np.int)
    reference_hist=np.array([0 for n in range(n_pinholes)], dtype=np.int)

    for k in range(attempts):
        randomized_photon_intensities = I0 * np.random.uniform(0., 1., n_pinholes)
    
        signal=np.array([0. for n in range(n_pinholes)], dtype=np.float)
        # E_S : signal photon: interacts with the pinhole detector
        reference=np.array([0. for n in range(n_pinholes)], dtype=np.float)
        # E_R : reference photon: interacts with the amplitude transmission mask
    
        idler_photon_directions = np.random.uniform(0., 180., n_pinholes)
        # Randomized idler photon direction
        for theta in idler_photon_directions:
            n = pinhole_number(theta)
            signal[n] += randomized_photon_intensities[n]
            reference[bucket_n[n]] += randomized_photon_intensities[n]
    
            signal_hist[n] = signal_hist[n] + 1
            reference_hist[bucket_n[n]] = reference_hist[bucket_n[n]] + 1
    
            bucket = np.sum(np.multiply(T, reference))
            # Interaction of photons and transmission mask
    
            correlated = bucket * np.array(signal)
            # Multiplying the bucket current and the original signal current
    
            correlated_signal += correlated
            # Sum up correlated signals
    
    plt.figure(i+1)
    plt.plot(correlated_signal/attempts, "ko")
    plt.ylabel("Ghost image [Average correlated intensities]")
    plt.xlabel("Transmission mask 1's degrees range [%.2f, %.2f]"%T_mask_range)

    print("Correlated signal: %s"%correlated_signal)
    print("Signal or idler photon intensities: %s"%signal)
    print("Reference photon intensities: %s"%reference)