# Plot intensity data

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
import pandas as pd
from qutip import *
import matplotlib.colors as colors
import matplotlib.cbook as cbook

from scipy.constants import *

from circuits import fluxonium
from scipy.constants import pi

rc('text', usetex=False)


### Long's earlier code
#Enter directory and name of measurement
# directory = 'G:\Projects\Fluxonium & qubits\Data\\2016_01\\13'
# measurement = 'S21_Phase2tones_ZNB_0&n40_BW100Hz_SMB_n20dBm_6to8GHz_YOKO_n7ton10V'
# path_data = directory + '\\' + measurement + '_Phase.csv'
# path_freq = directory + '\\' + measurement + '_Freq.csv'
# path_vol = directory + '\\' + measurement + '_Voltage.csv'
#
# RawData = np.genfromtxt(path_data, delimiter =',')
# Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
# V = np.genfromtxt(path_vol, delimiter =',')
#
# Z = RawData.transpose()
#
# #Optional: calculate differential
# Z_diff = np.diff(Z.transpose())
# Z_diff = Z_diff.transpose()
#
# #Plot the intensity map
# X, Y = np.meshgrid(V,Freq)
# fig=plt.figure()
# plt.pcolormesh(X, Y, Z, cmap=cm.BuGn, vmin =-1 , vmax =1)
# plt.title(measurement)
# plt.xlabel("Voltage (V)")
# plt.ylabel("Frequency (GHz)")
# plt.colorbar()
#
# #Click on the points on screen to define an approximation line for interpolation
# def onclick(event):
#     print '[%f, %f],'%(event.xdata, event.ydata)
# cid = fig.canvas.mpl_connect('button_press_event', onclick)
# plt.show()

##################### CHL code
# datadir_1 = 'Q3-current_sweep_pulsed_full_with_res_fitting_vs_flux_0.csv'
# datadir_1 = 'Q2-current_sweep_pulsed_full_with_res_fitting_vs_flux_5.csv'
datadir_1 = 'Q2-current_sweep_pulsed_full_with_res_fitting_vs_flux_fine.csv'
data = pd.read_csv(datadir_1, index_col=[0, 1])

# print('data=', data)


dfs = []
po_col = []
freq_range = []
Current = np.array([])

for C, sweep in data.groupby('Current'):
    # print(C)
    Current = np.append(Current, C)
    sweep = sweep.droplevel(0)

    logmag = sweep['amplitude']
    phase = sweep['phase']
    po_col += [[logmag, phase]]
    freq_range += [sweep.index]

Current = Current * 1e3
fig = plt.figure(figsize=(8, 6))
logmag_col = []
phase_col = []

logmag_col_normalized = []
phase_col_normalized = []

for i in range(np.shape(po_col)[0]):
    logmag_col += [po_col[i][0]]
    phase_col += [po_col[i][1]]
    logmag_col_normalized += [po_col[i][0] - np.median(po_col[i][0])]
    phase_col_normalized += [po_col[i][1] - np.median(po_col[i][1])]

Freq = np.array(freq_range[0]) / 1e9


X, Y = np.meshgrid(Current, Freq)
Z = np.array(np.transpose(logmag_col_normalized))


# plt.pcolormesh(X, Y, Z,  vmin=-1., vmax=1., cmap='RdBu_r')
# plt.pcolormesh(X, Y, Z,  vmin=-5., vmax=5.)
# plt.pcolormesh(X, Y, Z, vmin=-10., vmax=10.,
#                cmap='PuBu_r', shading='auto')

# plt.pcolormesh(X, Y, Z, vmin=-10., vmax=10.,
#                cmap='bwr', shading='auto')

plt.pcolormesh(X, Y, Z, vmin=-1., vmax=1.,
               cmap='bwr', shading='auto')

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
# plt.xlim([1, 8.5])
# plt.ylim([0, 7.2])
plt.xlabel("Current (mA)", fontsize=18)
plt.ylabel("Frequency (GHz)", fontsize=18)
# plt.colorbar()

# Click on the points on screen to define an approximation line for interpolation
def onclick(event):
    print('[%f, %f],' % (event.xdata, event.ydata))

cid = fig.canvas.mpl_connect('button_press_event', onclick)

if 0:   #(plot fitting data)
    E_J = 3.923616
    E_C = 0.968268
    E_L = 2.4814
    nlev_lc = 50
    nlev = 50

    I_o = (8.019 - 1.098) * 1e-3  # One period
    I_min = 4.565e-3  # half flux point current bias
    offset = (I_min - I_o / 2) / I_o  # 0 flux point

    phi_o = h / (2 * e)
    flux1 = Current * phi_o / I_o
    phi_ext1 = (flux1 / phi_o - offset) * 2 * np.pi

    phi_ext_array = np.linspace(0, 1, 101) * np.pi * 2

    qubit = fluxonium.Fluxonium_qubit(E_J=E_J, E_C=E_C, E_L=E_L,
                                      nlev=nlev, nlev_lc=nlev_lc)
    freq_01 = np.zeros_like(phi_ext_array)
    freq_12 = np.zeros_like(phi_ext_array)
    freq_02 = np.zeros_like(phi_ext_array)
    freq_03 = np.zeros_like(phi_ext_array)
    freq_13 = np.zeros_like(phi_ext_array)

    for phi_idx, phi_ext in enumerate(phi_ext_array):
        qubit.phi_ext = phi_ext
        freq_01[phi_idx] = qubit.freq(level1=0, level2=1)
        freq_12[phi_idx] = qubit.freq(level1=1, level2=2)
        freq_02[phi_idx] = qubit.freq(level1=0, level2=2)
        # freq_03[phi_idx] = qubit.freq(level1=0, level2=3)
        # freq_13[phi_idx] = qubit.freq(level1=1, level2=3)

    plt.plot(phi_ext_array / (2 * pi), freq_01)

if 0:
    E_j = 3.923616
    E_c = 0.968268
    E_l = 2.4814
    current = np.linspace(0, 10, 201) * 1e-3
    energy = np.zeros((len(current), 10))

    I_o = (8.019 - 1.098) * 1e-3  # One period
    I_min = 4.565e-3  # half flux point current bias
    offset = (I_min - I_o / 2) / I_o  # 0 flux point

    phi_o = h / (2 * e)
    flux1 = Current * phi_o / I_o
    phi_ext1 = (flux1 / phi_o - offset) * 2 * np.pi

    flux = current * phi_o / I_o
    phi_ext = (flux / phi_o - offset) * 2 * np.pi
    N = 50
    a = tensor(destroy(N))
    phi = (a + a.dag()) * (8.0 * E_c / E_l) ** 0.25 / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** 0.25 / np.sqrt(2.0)
    for idx in range(len(current)):
        ope = 1.0j * (phi - phi_ext[idx])
        H = 4.0 * E_c * na ** 2 + 0.5 * E_l * phi ** 2 - 0.5 * E_j * (ope.expm() + (-ope).expm())
        energy[idx, 0] = H.eigenenergies()[1] - H.eigenenergies()[0]
        energy[idx, 1] = H.eigenenergies()[2] - H.eigenenergies()[0]
        energy[idx, 2] = H.eigenenergies()[3] - H.eigenenergies()[0]
        energy[idx, 3] = H.eigenenergies()[2] - H.eigenenergies()[1]
        energy[idx, 4] = (H.eigenenergies()[2] - H.eigenenergies()[0])/2

    cut = 400
    plt.plot(current * 1e3, energy[:, 0], 'k--', dashes=(5, 5), label='01')
    plt.plot(current * 1e3, energy[:, 1], 'y--', dashes=(5, 5), label='02')
    # plt.plot(current * 1e3, energy[:, 2], '--', label='03')
    plt.plot(current * 1e3, energy[:, 3], 'r--', dashes=(5, 5), label='12')
    plt.plot(current * 1e3, energy[:, 4], 'g--', dashes=(5, 5), label='02/2')

    # plt.legend()

if 1:   #
    # data 'Q2-current_sweep_pulsed_full_with_res_fitting_vs_flux_5.csv'


    # E_l = 1.3007335909332736
    E_l = 1.26
    E_c = 0.9800357243876854
    E_j = 3.6834774239905403



    # current = np.linspace(-1.1, -0.7, 101) * 1e-3
    current = np.linspace(-1.6, -0.2, 101) * 1e-3
    energy = np.zeros((len(current), 10))

    I_o = 1.25 * 1e-3  # One period
    I_min = -0.89e-3  # half flux point current bias
    offset = (I_min - I_o / 2) / I_o  # 0 flux point

    phi_o = h / (2 * e)
    flux1 = Current * phi_o / I_o
    phi_ext1 = (flux1 / phi_o - offset) * 2 * np.pi

    flux = current * phi_o / I_o
    phi_ext = (flux / phi_o - offset) * 2 * np.pi
    N = 50
    a = tensor(destroy(N))
    phi = (a + a.dag()) * (8.0 * E_c / E_l) ** 0.25 / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** 0.25 / np.sqrt(2.0)
    for idx in range(len(current)):
        ope = 1.0j * (phi - phi_ext[idx])
        H = 4.0 * E_c * na ** 2 + 0.5 * E_l * phi ** 2 - 0.5 * E_j * (ope.expm() + (-ope).expm())
        energy[idx, 0] = H.eigenenergies()[1] - H.eigenenergies()[0]
        energy[idx, 1] = H.eigenenergies()[2] - H.eigenenergies()[0]
        energy[idx, 2] = H.eigenenergies()[3] - H.eigenenergies()[0]
        energy[idx, 3] = H.eigenenergies()[2] - H.eigenenergies()[1]
        energy[idx, 4] = (H.eigenenergies()[2] - H.eigenenergies()[0])/2

    cut = 400
    plt.plot(current * 1e3, energy[:, 0], 'k--', dashes=(5, 5), label='01')
    plt.plot(current * 1e3, energy[:, 1], 'y--', dashes=(5, 5), label='02')
    # plt.plot(current * 1e3, energy[:, 2], '--', label='03')
    plt.plot(current * 1e3, energy[:, 3], 'r--', dashes=(5, 5), label='12')
    plt.plot(current * 1e3, energy[:, 4], 'g--', dashes=(5, 5), label='02/2')

    E_l = 1.26
    E_c = 0.9800357243876854
    E_j = 3.6834774239905403


plt.title('Q2 Ej = 1.26 GHz, Ec = 0.98 GHz, El = 3.68 GHz (zoomed in)', fontsize=18)
plt.show()
