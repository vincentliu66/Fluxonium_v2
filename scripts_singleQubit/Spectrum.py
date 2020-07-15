"""
Plotting qubit spectrum

Author: Long Nguyen
Date: July 1, 2020
"""

import numpy as np
from matplotlib import pyplot as plt
from Fluxonium_v2.circuits import fluxonium

#Parameters
E_J = 3
E_C = 1
E_L = 1
nhilbert = 30
nlevels = 20
phi_ext_array = np.linspace(0,1,101)*np.pi*2
freqs = np.zeros((len(phi_ext_array), nlevels-1))
qubit = fluxonium.Fluxonium_qubit(E_J, E_C, E_L, 0, nhilbert, nlevels)
for phi_idx, phi_ext in enumerate(phi_ext_array):
    qubit.phi_ext = phi_ext
    for level_idx in range(1,nlevels):
        freqs[phi_idx, level_idx-1] = qubit.level(level_idx) - qubit.level(0)

for idx in range(nlevels - 1):
    plt.plot(phi_ext_array, freqs[idx])
plt.show()

