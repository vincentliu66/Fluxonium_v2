# -*- coding: utf-8 -*-
"""
Time evolution for fixed other parameters.

@author: Long Nguyen, Konstantin Nesterov
"""

import sys
sys.dont_write_bytecode = True
import time
import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate

import qutip as qt

import circuits.fluxonium as fluxonium
import circuits.coupled as coupled
import circuits.evolgates as evol
import devices
import scripts_singleQubit.plotting_settings

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


# Gate parameters.
T_gate = 100
transition_to_drive = ('11', '21')
# Scaling of the ideal value given by the inverse matrix element.
drive_amplitude_factor = 1  # 0.95436
# Drive frequency with respect to the resonance.
delta_omega_d = 0.0

# Pulse shape.
shape = 'gauss'  # 'gauss', 'cos' for 1-cos, or 'square'
sigma = 0.25  # sigma in units of T_gate for shape=='gauss'

# Method to calculate the propagator.
# 'propagator - full propagator using qt.propagator
# 'sesolve' - propagator using qt.sesolve for 4 computational states
method = 'propagator'

# Hilbert space.
# nlev_cav = 4
nlev_q = 6

save_figure = False
filename_prefix = 'stuff'

# Indices of the computational space.
comp_space = ['00', '01', '10', '11']
interaction = 'on'

# Transitions to show
states011 = ['11', '21', '22']
states010 = ['10', '20']
states001 = ['01', '02']
states000 = ['00']
###################################################################

time_start = time.time()

qubitA = fluxonium.Fluxonium_qubit(
    E_L=E_L1, E_C=E_C1, E_J=E_J1, nlev=nlev_q)
qubitB = fluxonium.Fluxonium_qubit(
    E_L=E_L2, E_C=E_C2, E_J=E_J2, nlev=nlev_q)
# resonator = cqed.Cavity(omega=inputs.omega_c, nlev=nlev_cav)

system = coupled.CoupledObjects(
    qubitA, qubitB, [qubitA, qubitB, J_C, 'charge'])
#            resonator, qubitA, qubitB,
#            [resonator, qubitA, inputs.g1, 'JC-charge'],
#            [resonator, qubitB, inputs.g2, 'JC-charge'],
#            [qubitA, qubitB, inputs.J_C, 'charge'])

level1, level2 = transition_to_drive[0], transition_to_drive[1]
# Calculate the drive frequency.
omega_d = abs(system.freq(level1, level2)) + delta_omega_d
# Calculate the drive amplitude.
matr_el = np.abs(system.n_ij(qubitA, level1, level2)
                 + system.n_ij(qubitB, level1, level2))
epsilon = drive_amplitude_factor / abs(matr_el)

print('Detuning between 10-20 and 11-21: {:.1f} MHz'.format(
    1000 * np.abs(system.freq('10', '20') - system.freq('11', '21'))))

print('Transition: {} - {} with frequency {:.4f} GHz'.format(
    level1, level2, abs(system.freq(level1, level2))))
print('Drive frequency: {:.4f} GHz'.format(omega_d))
print('Drive amplitude scale factor: {:.4f}'.format(
    drive_amplitude_factor))

t_points = np.linspace(0, T_gate, 2 * int(T_gate) + 1)
# The time-independent operator part of the drive term.
# a = system.a(resonator)
H_drive = epsilon * (system.n(0) + system.n(1))
# This calculates the evolution operator that works for
# computational levels only.
if method == 'sesolve':
    U_t = evol.evolution_compspace_microwave(
        system, H_drive, comp_space=comp_space, t_points=t_points,
        T_gate=T_gate, shape=shape, sigma=sigma, omega_d=omega_d,
        interaction=interaction)
elif method == 'propagator':
    U_t = evol.evolution_operator_microwave(
        system, H_drive, comp_space=comp_space, t_points=t_points,
        T_gate=T_gate, shape=shape, sigma=sigma, omega_d=omega_d,
        interaction=interaction)

fidelity = evol.fidelity_cz_gate(
    system, U_t, comp_space=comp_space,
    interaction=interaction)

print('max fidelity during the simulations: ', np.max(fidelity))

print('\n** Final values **')
print('Fidelity: ', fidelity[-1])
print('Diagonal elements of the evolution operator '
      + '(amplitudes and phases with respect to E*t in units of pi)')
U_f = U_t[-1]
U_me = {}
for state in comp_space:
    vec = system.eigvec(state, interaction=interaction)
    U_me[state] = U_f.matrix_element(vec.dag(), vec)
for state in comp_space:
    print(state, np.abs(U_me[state]),
          (np.angle(U_me[state]
                    * np.exp(2j * np.pi * system.level(state) * T_gate))) / np.pi)
print('(phi_00 + phi_11 - phi_01 - phi_01)/pi = ',
      (np.angle(U_me['00']) + np.angle(U_me['11'])
       - np.angle(U_me['01']) - np.angle(U_me['10'])) / np.pi)
phase_accum = (np.angle(U_me['00']) + np.angle(U_me['11'])
               - np.angle(U_me['01']) - np.angle(U_me['10']))
phase_accum = phase_accum / np.pi

P_driven_transition = evol.prob_transition(
    system, U_t, transition_to_drive[0],
    transition_to_drive[1], interaction=interaction)
t_2nd_excited = scipy.integrate.trapz(P_driven_transition, t_points)
print('Time spent in the 2nd state for {} - {}: {:.1f} ns'.format(
    transition_to_drive[0], transition_to_drive[1],
    t_2nd_excited))

P011 = {}
P010 = {}
P001 = {}
P000 = {}

# axes = np.empty((2, 2), dtype=object)
# fig = plt.figure(figsize=(12, 12))
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

# axes[0, 0] = fig.add_axes([0.1, 0.7, 0.35, 0.2])
# axes[0, 1] = fig.add_axes([0.6, 0.7, 0.35, 0.2])
# axes[1, 0] = fig.add_axes([0.1, 0.4, 0.35, 0.2])
# axes[1, 1] = fig.add_axes([0.6, 0.4, 0.35, 0.2])


ax011 = axes[0, 0]
ax010 = axes[0, 1]
ax001 = axes[1, 0]
ax000 = axes[1, 1]

for state in states011:
    P011[state] = evol.prob_transition(system, U_t, '11', state,
                                       interaction=interaction)
    ax011.plot(t_points, P011[state], lw=2,
               label=r'$P(11\rightarrow{})$'.format(state))

for state in states010:
    P010[state] = evol.prob_transition(
        system, U_t, '10', state, interaction=interaction)
    ax010.plot(t_points, P010[state], lw=2,
               label=r'$P(10\rightarrow {})$'.format(state))

for state in states001:
    P001[state] = evol.prob_transition(
        system, U_t, '01', state, interaction=interaction)
    ax001.plot(t_points, P001[state], lw=2,
               label=r'$P(01\rightarrow {})$'.format(state))

for state in states000:
    P000[state] = evol.prob_transition(
        system, U_t, '00', state, interaction=interaction)
    ax000.plot(t_points, P000[state], lw=2,
               label=r'$P(00\rightarrow {})$'.format(state))

textfontsize = 18
# fig.text(0.5, 0.3, r'At $t = {}$ ns: '.format(int(t_points[-1])),
#         fontsize=textfontsize, ha='center')
# fig.text(0.5, 0.25,
#         r'$P(11\rightarrow 11) = {:.4f}$, '.format(P011['11'][-1])
#         + r'$P(10\rightarrow 10) = {:.4f}$, '.format(P010['10'][-1])
#         + r'$P(01\rightarrow 10) = {:.4f}$, '.format(P001['01'][-1])
#         + r'$P(00\rightarrow 00) = {:.4f}$'.format(P000['00'][-1]),
#         fontsize=textfontsize, ha='center')
ax011.text(0.98, 0.93,
           r'$P(11 \rightarrow 11) = {:.6f}$'.format(P011['11'][-1]),
           ha='right', va='top', transform=ax011.transAxes,
           fontsize=textfontsize)
ax010.text(0.98, 0.93,
           r'$P(10 \rightarrow 10) = {:.6f}$'.format(P010['10'][-1]),
           ha='right', va='top', transform=ax010.transAxes,
           fontsize=textfontsize)
ax001.text(0.98, 0.93,
           r'$P(01 \rightarrow 01) = {:.6f}$'.format(P001['01'][-1]),
           ha='right', va='top', transform=ax001.transAxes,
           fontsize=textfontsize)
ax000.text(0.98, 0.93,
           r'$P(00 \rightarrow 00) = {:.6f}$'.format(P000['00'][-1]),
           ha='right', va='top', transform=ax000.transAxes,
           fontsize=textfontsize)
# fig.text(0.65, 0.25, r'$P(10) = {:.4f}$'.format(
#          P010['10'][-1]), fontsize=textfontsize)
fig.text(0.5, 0.1,
         r'CZ gate phase accumulation: '
         + '$\phi_{00} + \phi_{11} - \phi_{10} - \phi_{01} = $'
         + r'${:.3f} \pi $'.format(phase_accum),
         fontsize=textfontsize, ha='center')
fig.text(0.5, 0.05,
         r'Fidelity: '
         + r'$F = {:.6f}$'.format(fidelity[-1]),
         fontsize=textfontsize, ha='center')

for axarr in axes:
    for ax in axarr:
        ax.legend(loc='lower left')
        ax.set_xlim([np.min(t_points), np.max(t_points)])
        ax.set_xlabel('Time (ns)')
        ax.set_ylim([0, 1.02])
        ax.set_ylabel(r'$P_{i\rightarrow f}$')

ax011.set_title(
    r'Starting in $|11\rangle$')
ax010.set_title(
    r'Starting in $|10\rangle$')

ax001.set_title(
    r'Starting in $|01\rangle$')
ax000.set_title(
    r'Starting in $|00\rangle$')

fig.tight_layout(rect=[0, 0.15, 1, 1])


print('\nTotal time elapsed: ', time.time() - time_start)
plt.show()