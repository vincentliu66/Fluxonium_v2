"""
Plotting qubit spectrum

Author: Long Nguyen
Date: July 1, 2020
"""

import numpy as np
from matplotlib import pyplot as plt
from circuits import fluxonium
import plotting_settings


#Parameters
E_J = 3
E_C = 1
E_L = 1
nlev_lc = 30
nlev = 20
phi_ext_array = np.linspace(0,1,101)*np.pi*2

qubit = fluxonium.Fluxonium_qubit(E_J=E_J, E_C=E_C, E_L=E_L,
                                  nlev=nlev, nlev_lc=nlev_lc)

freq_01 = np.zeros_like(phi_ext_array)
freq_12 = np.zeros_like(phi_ext_array)
freq_02 = np.zeros_like(phi_ext_array)

for phi_idx, phi_ext in enumerate(phi_ext_array):
    qubit.phi_ext = phi_ext
    freq_01[phi_idx] = qubit.freq(level1 = 0, level2 = 1)
    freq_12[phi_idx] = qubit.freq(level1 = 1, level2 = 2)
    freq_02[phi_idx] = qubit.freq(level1 = 0, level2 = 2)

plt.plot(phi_ext_array, freq_01)
plt.plot(phi_ext_array, freq_12)
plt.plot(phi_ext_array, freq_02)
plt.show()

