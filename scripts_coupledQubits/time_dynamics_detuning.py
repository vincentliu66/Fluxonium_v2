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
from scipy.optimize import curve_fit

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
T_gate = 200
T_edge = 80
DRAG = True
DRAG_coefficient = 1.9
transition_to_drive = ('11', '21')
# Scaling of the ideal value given by the inverse matrix element.
drive_amplitude_factor = 0.0664*6.28  # 0.95436
# Drive frequency with respect to the resonance.
# delta_omega_d = 0.05
delta_omega_d_array = np.linspace(-0.0,0.0,1)

# Pulse shape.
shape = 'square'  # 'gauss', 'cos' for 1-cos, or 'square' or 'gauss_flat'
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
states011 = ['11', '21']
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
# '''
def oscillating_func(t,amp,freq,offset):
    return amp*np.cos(2*np.pi*freq*t)+offset

eff_freq = np.zeros_like(delta_omega_d_array)
eff_amp = np.zeros_like(delta_omega_d_array)
for wd_idx, delta_omega_d in enumerate(delta_omega_d_array):
    # Calculate the drive frequency.
    omega_d = abs(system.freq(level1, level2)) + delta_omega_d
    # Calculate the drive amplitude.
    matr_el = np.abs(system.n_ij(qubitA, level1, level2)
                     + system.n_ij(qubitB, level1, level2))
    epsilon = drive_amplitude_factor / abs(matr_el)

    t_points = np.linspace(0, T_gate, int(2*T_gate) + 1)
    # t_points = np.linspace(0, T_gate, 101)
    # t_points = np.linspace(0, 200, 2 * int(200) + 1)
    # The time-independent operator part of the drive term.
    H_drive = epsilon * (system.n(0) + system.n(1))
    U_t = evol.evolution_operator_microwave_long(
        system, H_drive, comp_space=comp_space, t_points=t_points,
        T_gate=T_gate, T_edge = T_edge, shape=shape, sigma=sigma,
        DRAG = DRAG, DRAG_coefficient = DRAG_coefficient,
        omega_d=omega_d, interaction=interaction)
    U_f = U_t[-1]
    U_me = {}
    for state in comp_space:
        vec = system.eigvec(state, interaction=interaction)
        U_me[state] = U_f.matrix_element(vec.dag(), vec)
    for state in comp_space:
        print(state, np.abs(U_me[state]),
              (np.angle(U_me[state]
                        * np.exp(2j * np.pi * system.level(state) * T_gate))) / np.pi)

    phase_accum = (np.angle(U_me['00']) + np.angle(U_me['11'])
                   - np.angle(U_me['01']) - np.angle(U_me['10']))
    phase_accum = phase_accum / np.pi

    P_driven_transition = evol.prob_transition(
        system, U_t, transition_to_drive[0],
        transition_to_drive[1], interaction=interaction)

    P011 = {}
    P010 = {}
    P001 = {}
    P000 = {}

    #Plotting
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    ax011 = axes[0]
    ax010 = axes[1]

    for state in states011:
        P011[state] = evol.prob_transition(system, U_t, '11', state,
                                           interaction=interaction)
        ax011.plot(t_points, P011[state], lw=2,
                   label=r'$P(11\rightarrow{})$'.format(state))

    # guess_amp = np.min(P011['21']) - np.max(P011['21'])
    # guess_freq = 0.06
    # guess = ([guess_amp, guess_freq,-1])
    # opt,cov = curve_fit(oscillating_func,xdata = t_points,ydata = P011['21'], p0=guess)
    eff_amp[wd_idx] = np.max(P011['21']) - np.min(P011['21'])
    # eff_freq[wd_idx] = opt[1]
    # ax011.plot(t_points, oscillating_func(t_points,*opt))

    for state in states010:
        P010[state] = evol.prob_transition(
            system, U_t, '10', state, interaction=interaction)
        ax010.plot(t_points, P010[state], lw=2,
                   label=r'$P_(10\rightarrow {})$'.format(state))

    print(str((wd_idx+1)/len(delta_omega_d_array)*100)+'%')
# '''
fname = '/Users/longnguyen/Documents/tmp'
# np.savetxt(fname+'_amp.txt',eff_amp)
# np.savetxt(fname+'_freq.txt',eff_freq)
textfontsize = 18

# eff_amp = np.genfromtxt(fname+'_amp.txt')
# eff_freq = np.genfromtxt(fname + '_freq.txt')

# ax011.text(0.98, 0.93,
#            r'$P(11 \rightarrow 11) = {:.6f}$'.format(P011['11'][-1]),
#            ha='right', va='top', transform=ax011.transAxes,
#            fontsize=textfontsize)
# ax010.text(0.98, 0.93,
#            r'$P(10 \rightarrow 10) = {:.6f}$'.format(P010['10'][-1]),
#            ha='right', va='top', transform=ax010.transAxes,
#            fontsize=textfontsize)

# for ax in axes:
#     ax.legend(loc='center right')
#     ax.set_xlim([np.min(t_points), np.max(t_points)])
#     ax.set_xlabel('Time (ns)')
#     ax.set_ylim([0, 1.02])
#     ax.set_ylabel(r'$P_{i\rightarrow f}$')
#
# ax011.set_title(
#     r'Starting in $|11\rangle$')
# ax010.set_title(
#     r'Starting in $|10\rangle$')

# fig.tight_layout(rect=[0, 0.15, 1, 1])

# detuning_fig, detuning_axes = plt.subplots(1, 2, figsize=(16, 9))
# # detuning_axes[0].plot(delta_omega_d_array*1e3, eff_freq*1e3 , 'd')
# # detuning_axes[1].plot(delta_omega_d_array*1e3, eff_amp, 'h')
# # detuning_axes[0].set_xlabel('Detuning (MHz)')
# # detuning_axes[1].set_xlabel('Detuning (MHz)')
# # detuning_axes[0].set_ylabel('Rabi freq (MHz)')
# # detuning_axes[1].set_ylabel('Rabi amplitude')

# detuning_fig = plt.figure(figsize=(16, 9))
# plt.plot(delta_omega_d_array*1e3, eff_amp, 'h')
# plt.xlabel('Detuning (MHz)')
# plt.ylabel('Leakage to 12')

plt.show()