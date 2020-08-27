#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stark-shift stuff vs omega_d.
"""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/konstantin/Work/work-sync/python_modules')
sys.path.insert(0, '../')

import numpy as np
from matplotlib import pyplot as plt
import time

import qutip as qt
import scripts_singleQubit.plotting_settings
import circuits.fluxonium as fluxonium
import circuits.coupled as coupled
import circuits.evolgates as evol
import devices


omega_d_range = np.linspace(4.49, 4.6, 11)

plt.close('all')

# Device parameters.
device_name = 'Augustus 17'
device = devices.devices[device_name]
E_L1 = device['parameters']['E_L1']
E_C1 = device['parameters']['E_C1']
E_J1 = device['parameters']['E_J1']
phi_ext1 = np.pi
E_L2 = device['parameters']['E_L2']
E_C2 = device['parameters']['E_C2']
E_J2 = device['parameters']['E_J2']
phi_ext2 = np.pi
coupling_type = device['coupling_type']
J_C = device['parameters']['J_C']

# Parameters for numerical calculations.
nlev_q_LC = 50  # The number of levels in the qubit LC-basis.
nlev_q = 5  # The number of qubit eigenstates.



##########################################################################
qubitA = fluxonium.Fluxonium_qubit(E_L=E_L1, E_C=E_C1, E_J=E_J1,
                        phi_ext=np.pi, nlev=nlev_q, nlev_lc=nlev_q_LC)
qubitB = fluxonium.Fluxonium_qubit(E_L=E_L2, E_C=E_C2, E_J=E_J2,
                        phi_ext=np.pi, nlev=nlev_q, nlev_lc=nlev_q_LC)

quant_sys = coupled.CoupledObjects(
        qubitA, qubitB, [qubitA, qubitB, J_C, coupling_type])

# ZZ coupling rate
xi_zz = quant_sys.freq('00', '10') - quant_sys.freq('01', '11')
print('ZZ coupling rate: {:.3f} MHz'.format(1000 * xi_zz))

T_stark = {}
for state in ['00', '01', '10', '11']:
    T_stark[state] = np.zeros_like(omega_d_range)
min_detuning = np.zeros_like(omega_d_range)
eps_th = np.zeros_like(omega_d_range)
T_total_stark = np.zeros_like(omega_d_range)

for ind, omega_d in enumerate(omega_d_range):
    print('ind, omega_d = ', ind, omega_d)
    min_detuning[ind] = 1000
    for i in range(5):
        for j in range(5):
            final = str(i) + str(j)
            for initial in ['00', '01', '10', '11']:
                matr_el2 = np.abs(quant_sys.n_ij(qubitA, initial, final)
                           + quant_sys.n_ij(qubitB, initial, final))**2
                freq_detuning = quant_sys.freq(initial, final) - omega_d
                T_stark[initial][ind] += -matr_el2 / freq_detuning

                if (matr_el2 > 1e-3 and
                    np.abs(freq_detuning) < min_detuning[ind]):
                    min_detuning[ind] = np.abs(freq_detuning)

T_total_stark = T_stark['01'] + T_stark['10'] - T_stark['00'] - T_stark['11']
for ind in range(len(omega_d_range)):
    if T_total_stark[ind] * xi_zz < 0:
        eps_th[ind] = 2 * np.sqrt(-xi_zz / T_total_stark[ind])

fig, axes = plt.subplots(2, figsize=(12, 12))

for state in ['00', '01', '10', '11']:
    axes[0].plot(omega_d_range, T_stark[state], label=state)

axes[1].plot(omega_d_range, T_total_stark)

axes[0].legend()

# for transition in [('00', '03'), ('00', '30'), ('10', '13'), ('01', '31'),
#                    ('01', '02'), ('10', '20'), ('11', '12'), ('11', '21')]:
#     print('{} - {} frequency: {:.3f} GHz'.format(transition[0], transition[1],
#             quant_sys.freq(transition[0], transition[1])))
#     axes[0].axvline(x=quant_sys.freq(transition[0], transition[1]),
#                     ls='--')
#     axes[1].axvline(x=quant_sys.freq(transition[0], transition[1]),
#                     ls='--')

for ax in [axes[0], axes[1]]:
    ax.set_ylim([-2, 2])

# if False:
#     axes[2].plot(omega_d_range, eps_th*1000, label='theoretical amplitude')
#     axes[2].plot(omega_d_range, min_detuning*1000, label='min detuning')
#     axes[2].legend()
#     axes[2].set_ylabel('MHz')
#
#     axes[3].plot(omega_d_range, eps_th / min_detuning, label='small parameter')
#     axes[3].legend()

for ax in axes:
    ax.set_xlim([np.min(omega_d_range), np.max(omega_d_range)])
    ax.set_xlabel(r'$\omega_d/2\pi$ (GHz)')

fig.tight_layout()
plt.show()
