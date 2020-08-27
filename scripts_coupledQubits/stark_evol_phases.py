#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AC Stark shift contribution: phases for ZZ rate.

"""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/konstantin/Work/work-sync/python_modules')
sys.path.insert(0, '../')

import numpy as np
from matplotlib import pyplot as plt
import time
import scipy.special


import qutip as qt

import cirqed as cqed
import myplotting
import devices
import evolgates_new as evol_new
import evolgates_old as evol


def main():
    #################### Define system parameters here #######################

    #device_name = 'AugustusVI'
    #device_name = 'Theory_Paper_2018'
    device_name = 'AugustusXVII'
    #device_name = 'AugustusXVIII'
    
    tmax = 100  # Maximum time.
    f_values = [0, 1, 2]#2.1, 2.15, 2.2]  # Drive amplitude values in units of eps_theoretical.
    
    delta_omega_d = 1  # Omega_d with respect to omega_00_21
    
    # Pulse shape.
    shape = 'two_phot_sin'
    sigma = 0.25 # sigma in units of T_gate for shape=='gauss'
    
    save_figure = False
    filename = device_name + '_negative_square' + '.png'
    
    # Parameters for numerical calculations.
    nlev_q_LC = 20  # The number of levels in the qubit LC-basis.                          
    nlev_q = 5  # The number of qubit eigenstates.
    
    ##########################################################################
    t0 = time.time()
    inputs = getattr(devices, device_name)()
        
    qubitA = cqed.Fluxonium(
            inputs.E_L1, inputs.E_C1, inputs.E_J1, nlev=nlev_q,
            nlev_lc=nlev_q_LC)
    qubitB = cqed.Fluxonium(
            inputs.E_L2, inputs.E_C2, inputs.E_J2, nlev=nlev_q,
            nlev_lc=nlev_q_LC)
    
    quant_system = cqed.CoupledObjects(
            qubitA, qubitB, [qubitA, qubitB, inputs.J_C, inputs.coupling])
    
    
    print('****************************************')
    print('Qubit A: \n')
    print('0-1: {:.3f} GHz'.format(qubitA.freq(0, 1)))
    print('1-2: {:.3f} GHz'.format(qubitA.freq(1, 2)))
    print('2-3: {:.3f} GHz'.format(qubitA.freq(2, 3)))
    print('0-2: {:.3f} GHz'.format(qubitA.freq(0, 2)))
    print('0-3: {:.3f} GHz'.format(qubitA.freq(0, 3)))
    print('Charge matrix element for the 0-1 transition: {:.3f}\n'.format(
            np.abs(qubitA.n_ij(0, 1))))
    print('Charge matrix element for the 1-2 transition: {:.3f}\n'.format(
            np.abs(qubitA.n_ij(1, 2))))
    print('Charge matrix element for the 0-3 transition: {:.3f}\n'.format(
            np.abs(qubitA.n_ij(0, 3))))
    
    print('Qubit B: \n')
    print('0-1: {:.3f} GHz'.format(qubitB.freq(0, 1)))
    print('1-2: {:.3f} GHz'.format(qubitB.freq(1, 2)))
    print('2-3: {:.3f} GHz'.format(qubitB.freq(2, 3)))
    print('0-2: {:.3f} GHz'.format(qubitB.freq(0, 2)))
    print('0-3: {:.3f} GHz'.format(qubitB.freq(0, 3)))
    print('Charge matrix element for the 0-1 transition: {:.3f}\n'.format(
            np.abs(qubitB.n_ij(0, 1))))
    print('Charge matrix element for the 1-2 transition: {:.3f}\n'.format(
            np.abs(qubitB.n_ij(1, 2))))
    print('Charge matrix element for the 0-3 transition: {:.3f}\n'.format(
            np.abs(qubitA.n_ij(0, 3))))
    
    omega_d = quant_system.freq('11', '21') + delta_omega_d
    
    t_points = np.linspace(0, tmax, 2 * int(tmax) + 1)
    
    # ZZ coupling rate
    xi_zz = quant_system.freq('00', '10') - quant_system.freq('01', '11')
    print('ZZ coupling rate: {:.3f} MHz'.format(1000 * xi_zz))
    
    # Get eps_theoretical.
    T_stark = {'00' : 0, '01' : 0, '10' : 0, '11' : 0}
    
    for i in range(5):
        for j in range(5):
            final = str(i) + str(j)
            for initial in ['00', '01', '10', '11']:
                matr_el2 = np.abs(quant_system.n_ij(qubitA, initial, final)
                           + quant_system.n_ij(qubitB, initial, final))**2
                freq_detuning = quant_system.freq(initial, final) - omega_d
                T_stark[initial] += -matr_el2 / freq_detuning

    T_total_stark = T_stark['01']+T_stark['10']-T_stark['00']-T_stark['11']
    print('T_stark : {:.2f} ns'.format(T_total_stark))
    eps_th = 2 * np.sqrt(np.abs(xi_zz / T_total_stark))
    print('Theoretical drive amplitude for correction: {:.2f} MHz'.format(
            1000 * eps_th))
    xi_cz = 0.01 / np.abs(quant_system.n_ij(qubitA, '11', '21')
            + quant_system.n_ij(qubitB, '11', '21'))
    print('Theoretical drive amplitude for 100 ns CZ gate: {:.2f} MHz'.format(
            1000 * xi_cz))
    
    print('Total time elapsed: ', time.time() - t0)
    phi_combination = np.zeros((len(f_values), len(t_points)))
    prob_combination = np.zeros_like(phi_combination)
    vec00 = quant_system.eigvec('00')
    vec01 = quant_system.eigvec('01')
    vec10 = quant_system.eigvec('10')
    vec11 = quant_system.eigvec('11')

    for ind, f in enumerate(f_values):
        t1 = time.time()
        print('Calculations for f = ', f)
    
        H_drive = tmax * f * eps_th *(quant_system.n(0) + quant_system.n(1))
        
        U_t = evol_new.evolution_subspace_microwave(
                    quant_system.H(), H_drive,
                    subspace_states=[quant_system.eigvec('00'), 
                    quant_system.eigvec('01'), quant_system.eigvec('10'), 
                    quant_system.eigvec('11')], 
                    t_points=t_points, 
                    T_gate=tmax, shape=shape, sigma=sigma, omega_d=omega_d)
        
        for tind in range(len(t_points)):
            u00 = U_t[tind].matrix_element(vec00.dag(), vec00)
            u01 = U_t[tind].matrix_element(vec01.dag(), vec01)
            u10 = U_t[tind].matrix_element(vec10.dag(), vec10)
            u11 = U_t[tind].matrix_element(vec11.dag(), vec11)
        
            phi_combination[ind, tind] = np.angle(u11*u00 / (u10*u01)) / np.pi
            prob_combination[ind, tind] = 0.25 * (np.abs(u11)**2 +
                    np.abs(u00)**2 + np.abs(u01)**2 + np.abs(u10)**2)
            
        print('Total time for this value of f: ', time.time() - t1)
    
    myplotting.my_figure_settings()
    fig, axes = plt.subplots(2, figsize=(10, 12))    

    for ind, f in enumerate(f_values):
        axes[0].plot(t_points, phi_combination[ind, :], 
            label=r'$\varepsilon/\varepsilon_{\rm th} = $' + '{:.2f}'.format(f))
        axes[1].plot(t_points, prob_combination[ind, :])
    
    for ax in axes:
        ax.set_xlim([0, tmax])
        ax.set_xlabel('Time (ns)')
        
    axes[1].set_ylim([0.96, 1.01])
    
    axes[0].legend()
    axes[0].set_title(r'$(\varphi_{11} + \varphi_{00}$'
                      + r'$- \varphi_{01} - \varphi_{10})/\pi$')
    axes[1].set_title(r'$(P_{00\to 00} + P_{01 \to 01} + P_{10\to 10}$'
                      +r'$ + P_{11\to 11})/4$')
    
    fig.text(0.5, 0.96, device_name + r', $\omega_{11\to 21}/2\pi = $'
            + '{:.3f} GHz'.format(quant_system.freq('11', '21'))
            + r', $\xi_{\rm ZZ}/2\pi = $'
            + r'{:.3f} MHz'.format(xi_zz * 1000), ha='center', fontsize=20)
    fig.text(0.5, 0.92,
            r'$(\omega_d -\omega_{11\to 21})/2\pi = $'
            + '{:.3f} GHz'.format(delta_omega_d)
            + r', $\varepsilon_{\rm th}/2\pi = $'
            + r'{:.1f} MHz'.format(eps_th * 1000)
            + r', $\varepsilon_{\rm CZ - 100 ns}/2\pi = $'
            + r'{:.1f} MHz'.format(1000 * xi_cz), ha='center', fontsize=20)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    
    if save_figure:
        fig.savefig(filename)
    
if __name__ == '__main__':
    main()
