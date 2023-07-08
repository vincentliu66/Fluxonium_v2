"""
Plotting qubit spectrum

Author: Long Nguyen
Date: July 1, 2020
"""

import numpy as np
from matplotlib import pyplot as plt
from circuits import fluxonium
import plotting_settings
from scipy.constants import pi


#Parameters
# E_J = 3.395
# E_C = 0.479
# E_L = 0.132

E_J = 5.4
E_C = 1.2
E_L = 1.5
nlev_lc = 50
nlev = 50
phi_ext_array = np.linspace(0,1,101)*np.pi*2

qubit = fluxonium.Fluxonium_qubit(E_J=E_J, E_C=E_C, E_L=E_L,
                                  nlev=nlev, nlev_lc=nlev_lc)

freq_01 = np.zeros_like(phi_ext_array)
freq_12 = np.zeros_like(phi_ext_array)
freq_23 = np.zeros_like(phi_ext_array)
freq_34 = np.zeros_like(phi_ext_array)
freq_02 = np.zeros_like(phi_ext_array)
freq_02Over2 = np.zeros_like(phi_ext_array)
freq_24 = np.zeros_like(phi_ext_array)
freq_13 = np.zeros_like(phi_ext_array)

freq_03Over3 = np.zeros_like(phi_ext_array)
freq_04Over4 = np.zeros_like(phi_ext_array)

for phi_idx, phi_ext in enumerate(phi_ext_array):
    qubit.phi_ext = phi_ext
    freq_01[phi_idx] = qubit.freq(level1 = 0, level2 = 1)
    freq_12[phi_idx] = qubit.freq(level1 = 1, level2 = 2)
    # freq_23[phi_idx] = qubit.freq(level1 = 2, level2 = 3)
    # freq_34[phi_idx] = qubit.freq(level1 = 3, level2 = 4)
    freq_02[phi_idx] = qubit.freq(level1=0, level2=2)
    freq_02Over2[phi_idx] = qubit.freq(level1=0, level2=2)/2

    freq_03Over3[phi_idx] = qubit.freq(level1=0, level2=3)/3
    freq_04Over4[phi_idx] = qubit.freq(level1=0, level2=4)/4

plt.plot(phi_ext_array/(2*pi), freq_03Over3, 'o', label='03/3')
plt.plot(phi_ext_array/(2*pi), freq_04Over4, 's', label='04/4')

# plt.plot(phi_ext_array/(2*pi), freq_23, 'o', label='23')
# plt.plot(phi_ext_array/(2*pi), freq_34, 's', label='34')
# plt.plot(phi_ext_array/(2*pi), freq_01, 'd', label='01')


plt.plot(phi_ext_array/(2*pi), freq_01, '-', label='01')
# plt.plot(phi_ext_array/(2*pi), freq_12, '--', label='12')
# plt.plot(phi_ext_array/(2*pi), freq_02, '+', label='02')
# plt.plot(phi_ext_array/(2*pi), freq_02Over2, 'x', label='02/2')

print(plt.rcParams['axes.prop_cycle'].by_key()['color'])
plt.legend()
plt.show()

