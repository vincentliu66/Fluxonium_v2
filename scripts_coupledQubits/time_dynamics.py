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
device_name = 'Augustus 17_2020/08'
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
transition_to_drive = ('11', '21')
# Scaling of the ideal value given by the inverse matrix element.
drive_amplitude_factor = 0.0664*2*np.pi  # 0.95436
# Drive frequency with respect to the resonance.
delta_omega_d = 0.1


# Pulse shape.
shape = 'gauss_flat_haonan'  # 'gauss', 'cos' for 1-cos, or 'square'
width = 10
T_flat = 100
T_gate = 2*width + T_flat
DRAG = True
DRAG_coefficient = 1.9

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

print ('10-20 freq:'+str(abs(system.freq('10', '20'))))
print ('11-21 freq:'+str(abs(system.freq('11', '21'))))

print ('nA(00-10)='+str(system.n_ij(qubitA, '00', '10')))
print ('nB(00-10)='+str(system.n_ij(qubitB, '00', '10')))

print ('nA(00-01)='+str(system.n_ij(qubitA, '00', '01')))
print ('nB(00-01)='+str(system.n_ij(qubitB, '00', '01')))

print ('nA(01-02)='+str(system.n_ij(qubitA, '01', '02')))
print ('nB(01-02)='+str(system.n_ij(qubitB, '01', '02')))

print ('nA(11-12)='+str(system.n_ij(qubitA, '11', '12')))
print ('nB(11-12)='+str(system.n_ij(qubitB, '11', '12')))

print ('nA(10-20)='+str(system.n_ij(qubitA, '10', '20')))
print ('nB(10-20)='+str(system.n_ij(qubitB, '10', '20')))

print ('nA(11-21)='+str(system.n_ij(qubitA, '11', '21')))
print ('nB(11-21)='+str(system.n_ij(qubitB, '11', '21')))

#30 <-> 13
print ('nA(20-30)='+str(system.n_ij(qubitA, '20', '13')))
print ('nB(20-30)='+str(system.n_ij(qubitB, '20', '13')))

print ('nA(21-31)='+str(system.n_ij(qubitA, '21', '31')))
print ('nB(21-31)='+str(system.n_ij(qubitB, '21', '31')))

print ('nA(00-03)='+str(system.n_ij(qubitA, '00', '03')))
print ('nB(00-03)='+str(system.n_ij(qubitB, '00', '03')))

#30 <-> 13
print ('nA(10-13)='+str(system.n_ij(qubitA, '10', '30')))
print ('nB(10-013)='+str(system.n_ij(qubitB, '10', '30')))

print ('nA(01-31)='+str(system.n_ij(qubitA, '01', '31')))
print ('nB(01-31)='+str(system.n_ij(qubitB, '01', '31')))

#30 <-> 13
print ('nA(00-30)='+str(system.n_ij(qubitA, '00', '13')))
print ('nB(00-30)='+str(system.n_ij(qubitB, '00', '13')))

print ('nA(00-31)='+str(system.n_ij(qubitA, '00', '31')))
print ('nB(00-31)='+str(system.n_ij(qubitB, '00', '31')))

#30 <-> 13
print ('nA(00-13)='+str(system.n_ij(qubitA, '00', '30')))
print ('nB(00-13)='+str(system.n_ij(qubitB, '00', '30')))

# Calculate the drive frequency.
omega_d = abs(system.freq(level1, level2)) + delta_omega_d
# Calculate the drive amplitude.
matr_el = np.abs(system.n_ij(qubitA, level1, level2)
                 + system.n_ij(qubitB, level1, level2))
epsilon = drive_amplitude_factor / abs(matr_el)

t_points = np.linspace(0, T_gate, 2 * int(T_gate) + 1)
# t_points = np.linspace(0, T_gate, 101)
# t_points = np.linspace(0, 200, 2 * int(200) + 1)
# The time-independent operator part of the drive term.
H_drive = epsilon * (system.n(0) + system.n(1))
U_t = evol.evolution_operator_microwave_long(
    system, H_drive, comp_space=comp_space, t_points=t_points,
    T_flat=T_flat, T_gate=T_gate, shape=shape, width=width, omega_d=omega_d,
    DRAG = DRAG, DRAG_coefficient = DRAG_coefficient,
    interaction=interaction)
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

for state in states010:
    P010[state] = evol.prob_transition(
        system, U_t, '10', state, interaction=interaction)
    ax010.plot(t_points, P010[state], lw=2,
               label=r'$P(10\rightarrow {})$'.format(state))

textfontsize = 18
ax011.text(0.98, 0.93,
           r'$P(11 \rightarrow 11) = {:.6f}$'.format(P011['11'][-1]),
           ha='right', va='top', transform=ax011.transAxes,
           fontsize=textfontsize)
ax010.text(0.98, 0.93,
           r'$P(10 \rightarrow 10) = {:.6f}$'.format(P010['10'][-1]),
           ha='right', va='top', transform=ax010.transAxes,
           fontsize=textfontsize)


for ax in axes:
    ax.legend(loc='lower left')
    ax.set_xlim([np.min(t_points), np.max(t_points)])
    ax.set_xlabel('Time (ns)')
    ax.set_ylim([0, 1.02])
    ax.set_ylabel(r'$P_{i\rightarrow f}$')

ax011.set_title(
    r'Starting in $|11\rangle$')
ax010.set_title(
    r'Starting in $|10\rangle$')

fig.tight_layout(rect=[0, 0.15, 1, 1])
plt.show()