"""
Created on 2020/08/11

@author: Long Nguyen, based on code by Konstantin Nesterov
"""

import sys

sys.dont_write_bytecode = True

import time
import numpy as np
from matplotlib import pyplot as plt
import scripts_singleQubit.plotting_settings
import devices
import circuits.fluxonium as fluxonium
import circuits.coupled as coupled

plt.close("all")
##############################################################################

time_start = time.time()

take_data_from_input_file = True
device_name = 'Augustus 17'
nlev_single = 7  # The number of single-qubit levels to show.
nlev_show = 15  # The number of two-qubit levels to show.

if not take_data_from_input_file:
    # Parameters of the first fluxonium.
    E_L1 = 0.83  # inductive energy
    E_C1 = 1.03  # charging energy
    E_J1 = 5.35  # Josephson energy
    phi_ext1 = np.pi  # external phase shift

    # Parameters of the second fluxonium.
    E_L2 = 0.81  # inductive energy
    E_C2 = 1.06  # charging energy
    E_J2 = 3.93  # Josephson energy
    phi_ext2 = np.pi  # external phase shift

    # Interaction energy between two fluxoniums.
    # E_int n_1 n_2 or E_int phi_1 phi_2.
    E_int = 0.27  # when a single value is needed
    E_int_range = np.linspace(0, 0.3, 11)  # when a range is needed
    coupling = 'charge'  # 'charge' or 'flux'
else:
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
    if coupling_type == 'charge':
        E_int = device['parameters']['J_C']
        E_int_range = np.linspace(0,E_int*1.2,11)
    else:
        E_int = device['parameters']['J_L']
        E_int_range = np.linspace(0,E_int*1.2,11)

qubit1 = fluxonium.Fluxonium_qubit(E_J = E_J1, E_C=E_C1, E_L=E_L1, phi_ext=phi_ext1)
qubit2 = fluxonium.Fluxonium_qubit(E_J = E_J2, E_C=E_C2, E_L=E_L2, phi_ext=phi_ext2)

#######################################################################################################################
#######################################################################################################################
# Get spectra of two bare fluxoniums as a function of phi_ext.
phi_points = np.linspace(0, 2 * np.pi, 100)

energies1 = np.zeros((nlev_single, len(phi_points)))
energies2 = np.zeros_like(energies1)
for iphi, phi_ext in enumerate(phi_points):
    qubit1.phi_ext = phi_ext
    qubit2.phi_ext = phi_ext
    energies1[:, iphi] = qubit1.levels(nlev=nlev_single)
    energies2[:, iphi] = qubit2.levels(nlev=nlev_single)

qubit1.phi_ext = phi_ext1
qubit2.phi_ext = phi_ext2

# Figure for the noninteracting qubits
fig_nonint = plt.figure(figsize=(10, 10), dpi=100)
fig_nonint.suptitle('Bare fluxoniums', fontsize=20)
axes_fl = []
axes_fl.append(fig_nonint.add_axes([0.1, 0.5, 0.35, 0.35]))
axes_fl.append(fig_nonint.add_axes([0.55, 0.5, 0.35, 0.35]))

for idx in range(7):
    axes_fl[0].plot(phi_points, energies1[idx, :], linewidth=2)
    axes_fl[1].plot(phi_points, energies2[idx, :], linewidth=2)

axes_fl[0].set_title('Qubit 1', fontsize=20)
axes_fl[1].set_title('Qubit 2', fontsize=20)

for idx in range(2):
    axes_fl[idx].set_xlim([-0.2, 6.4])
    axes_fl[idx].set_xticks([0, np.pi, 2 * np.pi])
    axes_fl[idx].set_xticklabels(['0', r'$\pi$', r'$2\pi$'], fontsize=18)
    axes_fl[idx].set_xlabel(r'$\varphi_{\rm ext}$', fontsize=18)
    axes_fl[idx].set_ylim([-2, 20.2])
    axes_fl[idx].set_ylabel('Energy (GHz)', fontsize=18)
    axes_fl[idx].tick_params(labelsize = 18.0)
    axes_fl[idx].grid(True)

string1 = (r'$E_L = $' + str(round(E_L1, 2)) + ' GHz, '
           + r'$E_C = $' + str(round(E_C1, 2)) + ' GHz, \n'
           + r'$E_J = $' + str(round(E_J1, 2)) + ' GHz \n\n'
           + r'At $\varphi_{\rm ext}/2\pi = $'
           + str(round(phi_ext1 / np.pi/2, 1)) + ':\n '
           + r'$\nu_{01} = $'
           + str(round(qubit1.freq(level1=0, level2=1), 2))
           + ' GHz,\n ' + r'$\nu_{12} = $'
           + str(round(qubit1.freq(level1=1, level2=2), 2))
           + ' GHz,\n ' + r'$\nu_{03} = $'
           + str(round(qubit1.freq(level1=0, level2=3), 2))
           + ' GHz')
string2 = (r'$E_L = $' + str(round(E_L2, 2)) + ' GHz, '
           + r'$E_C = $' + str(round(E_C2, 2)) + ' GHz, \n'
           + r'$E_J = $' + str(round(E_J2, 2)) + ' GHz \n\n'
           + r'At $\varphi_{\rm ext}/2\pi = $'
           + str(round(phi_ext2 / np.pi/2, 1)) + ':\n '
           + r'$\nu_{01} = $'
           + str(round(qubit2.freq(level1=0, level2=1), 2))
           + ' GHz,\n ' + r'$\nu_{12} = $'
           + str(round(qubit2.freq(level1=1, level2=2), 2))
           + ' GHz,\n ' + r'$\nu_{03} = $'
           + str(round(qubit2.freq(level1=0, level2=3), 2))
           + ' GHz')
fig_nonint.text(0.1, 0.2, string1, fontsize=18)
fig_nonint.text(0.55, 0.2, string2, fontsize=18)

#######################################################################################################################
#######################################################################################################################
#Interacting system
energies_coupled = np.zeros(
    (qubit1.nlev * qubit2.nlev, len(E_int_range)))
n1_11_21 = np.zeros(len(E_int_range), dtype=complex)
n2_11_21 = np.zeros_like(n1_11_21, dtype=complex)
n1_10_20 = np.zeros_like(n1_11_21, dtype=complex)
n2_10_20 = np.zeros_like(n1_11_21, dtype=complex)
for indE, E_int in enumerate(E_int_range):
    system = coupled.CoupledObjects(qubit1, qubit2,
                                    [qubit1, qubit2, E_int, coupling_type])
    energies_coupled[:, indE] = system.levels()
    n1_11_21[indE] = system.n_ij(qubit1, '11', '21')
    n1_10_20[indE] = system.n_ij(qubit1, '10', '20')
    n2_11_21[indE] = system.n_ij(qubit2, '11', '21')
    n2_10_20[indE] = system.n_ij(qubit2, '10', '20')

print('0-1 frequencies detuning (omega_A - omega_B): ',
      qubit1.freq(0, 1) - qubit2.freq(0, 1), ' GHz\n')
print('Noninteracting detuning between 21 and 12 (E_12 - E_21): ',
      system.freq_nonint('21', '12'),
      ' GHz\n')
print('Time elapsed: ', time.time() - time_start)

# Figure for the interacting qubits
fig_coupled = plt.figure(figsize=(16, 9), dpi=100)
fig_coupled.suptitle('Energies of two interacting qubits',
                     fontsize=24)
axes_coupled1 = fig_coupled.add_axes([0.05, 0.1, 0.25, 0.7])
axes_coupled2 = fig_coupled.add_axes([0.35, 0.1, 0.25, 0.7])
axes_coupled3 = fig_coupled.add_axes([0.65, 0.1, 0.25, 0.7])

if coupling_type == 'charge':
    string_temp = r'$E_{\rm int} \hat{n}_1 \hat{n}_2$ interaction'
else:
    string_temp = r'$E_{\rm int} \hat{\phi}_1 \hat{\phi}_2$ interaction'
fig_coupled.text(0.5, 0.9, string_temp, fontsize=20, ha='center')

string_temp = ('at ' + r'$\varphi_{\rm ext1}/\pi = $'
               + str(round(phi_ext1 / np.pi, 1))
               + ' and ' + r'$\varphi_{\rm ext2}/\pi = $'
               + str(round(phi_ext2 / np.pi, 1)))
fig_coupled.text(0.5, 0.85, string_temp, fontsize=20, ha='center')

for idx in range(nlev_show - 1, -1, -1):
    label = (r'$|$' + system.level_label(idx, label_format='str') + r'$\rangle$')
    axes_coupled1.plot(E_int_range, energies_coupled[idx, :],
                       linewidth=2, label=label)
    axes_coupled2.text(0, 0.01 + 0.05 * idx, label, fontsize=12)
    axes_coupled2.plot(E_int_range, energies_coupled[idx, :] \
                       - energies_coupled[idx, 0] + 0.05 * idx, linewidth=2)
    axes_coupled3.plot(E_int_range, energies_coupled[idx, :] - energies_coupled[0, :], linewidth =2, label = label)


axes_coupled1.set_xlabel(r'$E_{\rm int}$ (GHz)',
                         fontsize=18)
axes_coupled1.set_ylabel(r'$E$ (GHz)', fontsize=18)
axes_coupled1.grid(ls='--', lw=1)
axes_coupled1.legend(loc=2, prop={'size':12})
# axes_coupled1.tick_params(labelsize=16.0)
axes_coupled2.set_xlabel(r'$E_{\rm int}$ (GHz)',
                         fontsize=18)
axes_coupled2.set_ylabel(r'$E(E_{\rm int}) - E(0)$ (GHz)', fontsize=18)
axes_coupled2.grid(ls='--', lw=1)
axes_coupled2.set_title('Changes in energies (with vertical offsets)', size = 14.0)
# axes_coupled2.tick_params(labelsize=16.0)
axes_coupled3.set_title('Trans energies (GHz)', size = 14.0)
axes_coupled3.set_ylabel(r'$E-E_{00}$ (GHz)', fontsize=18)
axes_coupled3.legend(loc=2, prop={'size':12})
# axes_coupled3.tick_params(labelsize=16.0)
#######################################################################################################################
#######################################################################################################################
# Figure for energies in the computational space.

ind00 = system.level_label('00')
ind01 = system.level_label('01')
ind10 = system.level_label('10')
ind11 = system.level_label('11')
ind20 = system.level_label('20')
ind21 = system.level_label('21')

comb11 = (energies_coupled[ind11, :] + energies_coupled[ind00, :]
          - energies_coupled[ind01, :] - energies_coupled[ind10, :])
comb21 = (energies_coupled[ind21, :] + energies_coupled[ind00, :]
          - energies_coupled[ind01, :] - energies_coupled[ind20, :])
comb_new = (energies_coupled[ind21, :] - energies_coupled[ind11, :]
            - energies_coupled[ind20, :] + energies_coupled[ind10, :])

fig_comb3, axes_comb3 = plt.subplots(figsize=(16, 9), dpi=100)
axes_comb3.plot(E_int_range,
                1000 * (energies_coupled[ind00, :] - energies_coupled[ind00, 0]),
                lw=2, ls='--', label=r'$\Delta E_{|00\rangle}$')
axes_comb3.plot(E_int_range,
                1000 * (energies_coupled[ind10, :] - energies_coupled[ind10, 0]),
                lw=2, ls='-.', label=r'$\Delta E_{|10\rangle}$')
axes_comb3.plot(E_int_range,
                1000 * (energies_coupled[ind01, :] - energies_coupled[ind01, 0]),
                lw=2, ls='-.', label=r'$\Delta E_{|01\rangle}$')
axes_comb3.plot(E_int_range,
                1000 * (energies_coupled[ind11, :] - energies_coupled[ind11, 0]),
                lw=2, ls=':', label=r'$\Delta E_{|11\rangle}$')
axes_comb3.plot(E_int_range, 1000 * comb11, lw=3, ls='-',
                label=r'$E_{|11\rangle} + E_{|00\rangle}$'
                      + r'$- E_{|01\rangle} - E_{|10\rangle}$')
axes_comb3.set_xlabel(r'$E_{\rm int}$ (GHz)', fontsize=18)
axes_comb3.set_ylabel('Energy (MHz)', fontsize=18)
axes_comb3.grid(ls='--', lw=1)
axes_comb3.legend(fontsize=16)
# axes_comb3.set_ylim([-6, 0])
axes_comb3.set_xlim([min(E_int_range), max(E_int_range)])
axes_comb3.set_title('What happens in the computational space.')
axes_comb3.tick_params(labelsize = 20)
#######################################################################################################################
#######################################################################################################################
# Figure for the matrix elements.

fig_me = plt.figure(figsize=(16, 9), dpi=100)
axes_me = fig_me.add_axes([0.1, 0.1, 0.8, 0.7])

axes_me.plot(E_int_range, abs(n1_11_21), lw=2, ls='-',
             label=r'$|\langle 11 | \hat{n}_1 | 21\rangle| $')
axes_me.plot(E_int_range, abs(n1_10_20), lw=2, ls='--',
             label=r'$|\langle 10 | \hat{n}_1 | 20\rangle| $')
axes_me.plot(E_int_range, abs(n2_11_21), lw=2, ls='-.',
             label=r'$|\langle 11 | \hat{n}_2 | 21\rangle| $')
axes_me.plot(E_int_range, abs(n2_10_20), lw=2, ls=':',
             label=r'$|\langle 10 | \hat{n}_2 | 20\rangle| $')

axes_me.legend(fontsize=16)
axes_me.set_title('Matrix elements of the charge operators', fontsize=20)
axes_me.set_xlabel(r'$E_{\rm int}$ (GHz)', fontsize=18)
axes_me.grid(ls='--', lw=1)
axes_me.tick_params(labelsize = 20)
#######################################################################################################################
#######################################################################################################################
# Figure for the always-on interaction

fig_proposal_new, axes_proposal_new = plt.subplots(figsize=(16, 9), dpi=100)

axes_proposal_new.plot(1000 * E_int_range, 100000 * np.abs(comb11), lw=2, ls='-',
                       label=r'$|E_{|\widetilde{11}\rangle} - E_{|\widetilde{10}\rangle}$'
                             + r'$ - E_{|\widetilde{01}\rangle} + E_{|\widetilde{00}\rangle}|$')
axes_proposal_new.plot(1000 * E_int_range, 1000 * np.abs(comb_new), lw=2, ls='--',
                       label=r'$|E_{|\widetilde{21}\rangle} - E_{|\widetilde{11}\rangle}$'
                             + r'$- \left(E_{|\widetilde{20}\rangle} - E_{|\widetilde{10}\rangle}\right)|$')

axes_proposal_new.set_ylabel(r'Differential Energy (MHz)', fontsize=20)
axes_proposal_new.legend(loc=3, fontsize=16)

axes_proposal_new.set_xlabel(r'$E_{\rm int}/h$ (MHz)', fontsize=20)
axes_proposal_new.grid(ls='--', lw=1)
axes_proposal_new.set_xlim(left=0, right=1000 * max(E_int_range))
axes_proposal_new.tick_params(length=7, width=2, direction='in', right='on',
                              top='on', labelsize=16)
for axis in ['left', 'right', 'top', 'bottom']:
    axes_proposal_new.spines[axis].set_linewidth(2)
ind = int(0.9 * len(E_int_range))
axes_proposal_new.annotate(
    r'$\times 100$', xy=(1000 * E_int_range[ind], 100000 * np.abs(comb11[ind])),
    xytext=(-100, -10), textcoords='offset points',
    fontsize=20, ha='center', va='center',
    arrowprops=dict(arrowstyle='->', lw=2, mutation_scale=15))

########################################################################################################################
print('total time elapsed for `spectrum_double.py` = ',
          time.time() -time_start, 's\n')
plt.show()