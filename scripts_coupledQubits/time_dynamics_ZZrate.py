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
import circuits.evolgates_new
import devices
import scripts_singleQubit.plotting_settings

plt.close('all')

def oscillating_func(t,amp,freq,offset):
    return amp*np.cos(2*np.pi*freq*t)+offset

# Device parameters.
device_name = 'Augustus 17_2020/08'
device = devices.devices[device_name]
E_L1 = device['parameters']['E_L1']
E_C1 = device['parameters']['E_C1']
E_J1 = device['parameters']['E_J1']
phi_ext1 = 1.01*np.pi
E_L2 = device['parameters']['E_L2']
E_C2 = device['parameters']['E_C2']
E_J2 = device['parameters']['E_J2']
phi_ext2 = 1.01*np.pi
coupling_type = device['coupling_type']
J_C = device['parameters']['J_C']

# Scaling of the ideal value given by the inverse matrix element.
drive_amplitude_factor_array = np.linspace(0.2, 0.7, 51)*2*np.pi*0.0524/0.5
delta_omega_d_array = np.linspace(0.01,0.08,71)

# Pulse shape.
shape = 'gauss_flat_haonan'  # 'gauss', 'cos' for 1-cos, or 'square' or 'gauss_flat'
sigma = 0.25  # sigma in units of T_gate for shape=='gauss'
DRAG = False
DRAG_coefficient = 1.9
transition_to_drive = ('11', '21')
width = 20
T_flat = 160
T_gate = T_flat+2*width

# Method to calculate the propagator.
# 'propagator - full propagator using qt.propagator
# 'sesolve' - propagator using qt.sesolve for 4 computational states
method = 'sesolve'

# Hilbert space.
# nlev_cav = 4
nlev_q = 7

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
    E_L=E_L1, E_C=E_C1, E_J=E_J1, nlev=nlev_q, phi_ext=phi_ext1)
qubitB = fluxonium.Fluxonium_qubit(
    E_L=E_L2, E_C=E_C2, E_J=E_J2, nlev=nlev_q, phi_ext=phi_ext2)
# resonator = cqed.Cavity(omega=inputs.omega_c, nlev=nlev_cav)
system = coupled.CoupledObjects(
    qubitA, qubitB, [qubitA, qubitB, J_C, 'charge'])
#            resonator, qubitA, qubitB,
#            [resonator, qubitA, inputs.g1, 'JC-charge'],
#            [resonator, qubitB, inputs.g2, 'JC-charge'],
#            [qubitA, qubitB, inputs.J_C, 'charge'])

# print('00-03: '+str(abs(system.freq('00', '03'))))
# print('00-30: '+str(abs(system.freq('00', '30'))))
# print('00-13: '+str(abs(system.freq('00', '13'))))
# print('00-31: '+str(abs(system.freq('00', '31'))))
# print('01-31: '+str(abs(system.freq('01', '31'))))
# print('10-13: '+str(abs(system.freq('10', '13'))))
# print('11-21: '+str(abs(system.freq('11', '21'))))

level1, level2 = transition_to_drive[0], transition_to_drive[1]
coupl_ratio = 0.7718
t_points = np.linspace(0, T_gate, 2*int(T_gate) + 1)
phase = np.zeros((len(drive_amplitude_factor_array), len(delta_omega_d_array), len(t_points)))
for a_idx, drive_amplitude_factor in enumerate(drive_amplitude_factor_array):
    for dwd_idx, delta_omega_d in enumerate(delta_omega_d_array):
        # Calculate the drive frequency.
        omega_d = abs(system.freq(level1, level2)) + delta_omega_d
        # Calculate the drive amplitude.
        matr_el = np.abs(system.n_ij(qubitA, level1, level2)
                         + coupl_ratio*system.n_ij(qubitB, level1, level2))
        epsilon = drive_amplitude_factor / abs(matr_el)
        # The time-independent operator part of the drive term.
        H_drive = epsilon * (system.n(0) + coupl_ratio*system.n(1))
        if method == 'sesolve':
            U_t = evol.evolution_compspace_microwave_long(
                system, H_drive, comp_space=comp_space, t_points=t_points,
                T_gate=T_gate, T_flat=T_flat, width = width, shape=shape,
                DRAG = DRAG, omega_d=omega_d, interaction=interaction)
        elif method == 'propagator':
            U_t = evol.evolution_operator_microwave_long(
                system, H_drive, comp_space=comp_space, t_points=t_points,
                T_gate=T_gate, T_flat=T_flat, width = width, shape=shape,
                DRAG=DRAG, omega_d=omega_d, interaction=interaction)

        vec00 = system.eigvec('00')
        vec01 = system.eigvec('01')
        vec10 = system.eigvec('10')
        vec11 = system.eigvec('11')
        for tind in range(len(t_points)):
            u00 = U_t[tind].matrix_element(vec00.dag(), vec00)
            u01 = U_t[tind].matrix_element(vec01.dag(), vec01)
            u10 = U_t[tind].matrix_element(vec10.dag(), vec10)
            u11 = U_t[tind].matrix_element(vec11.dag(), vec11)

            phase[a_idx, dwd_idx, tind] = np.angle(u11 * u00 / (u10 * u01))

        completed = dwd_idx+1 + len(drive_amplitude_factor_array)*a_idx
        print(str(np.round(completed/ (len(drive_amplitude_factor_array)*len(drive_amplitude_factor_array)) * 100, 2)) + '% completed')
        print (str(np.round(time.time()-time_start, 3) )+ 's has elapsed')

fname = '/Users/longnguyen/Documents/tmp'
np.save(fname+'_phaseZZ_1121_away',phase)

# plt.figure(figsize = [16,9])
# plt.plot(t_points, np.unwrap(phase[0,0,:]),'-h')
# # print ('time elapsed = ' + str(time.time()-time_start))
# plt.xlabel('ns')
# plt.ylabel('phase (rad)')
# plt.show()