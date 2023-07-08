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


##################### CHL code
datadir_1 = 'Q4-zoomed_in_current_sweep_pulsed_full_with_res_fitting_vs_flux_2.csv'
datadir_2 = 'Q4-zoomed_in_current_sweep_pulsed_full_with_res_fitting_vs_flux_3.csv'
data1 = pd.read_csv(datadir_1, index_col=[0, 1])
data2 = pd.read_csv(datadir_2, index_col=[0, 1])

# print('data=', data)

if 1:
    dfs1 = []
    po_col1 = []
    freq_range1 = []
    Current1 = np.array([])

    for C, sweep in data1.groupby('Current'):
        # print(C)
        Current1 = np.append(Current1, C)
        sweep = sweep.droplevel(0)

        logmag1 = sweep['amplitude']
        phase1 = sweep['phase']
        po_col1 += [[logmag1, phase1]]
        freq_range1 += [sweep.index]

    Current1 = Current1 * 1e3
    fig = plt.figure(figsize=(8, 6))
    logmag_col1 = []
    phase_col1 = []

    logmag_col_normalized1 = []
    phase_col_normalized1 = []

    for i in range(np.shape(po_col1)[0]):
        logmag_col1 += [po_col1[i][0]]
        phase_col1 += [po_col1[i][1]]
        logmag_col_normalized1 += [po_col1[i][0] - np.median(po_col1[i][0])]
        phase_col_normalized1 += [po_col1[i][1] - np.median(po_col1[i][1])]

    Freq1 = np.array(freq_range1[0]) / 1e9


    X1, Y1 = np.meshgrid(Current1, Freq1)
    Z1 = np.array(np.transpose(logmag_col_normalized1))

if 1:
    dfs2 = []
    po_col2 = []
    freq_range2 = []
    Current2 = np.array([])

    for C, sweep in data2.groupby('Current'):
        # print(C)
        Current2 = np.append(Current2, C)
        sweep = sweep.droplevel(0)

        logmag2 = sweep['amplitude']
        phase2 = sweep['phase']
        po_col2 += [[logmag2, phase2]]
        freq_range2 += [sweep.index]

    Current2 = Current2 * 1e3
    fig = plt.figure(figsize=(8, 6))
    logmag_col2 = []
    phase_col2 = []

    logmag_col_normalized2 = []
    phase_col_normalized2 = []

    for i in range(np.shape(po_col2)[0]):
        logmag_col2 += [po_col2[i][0]]
        phase_col2 += [po_col2[i][1]]
        logmag_col_normalized2 += [po_col2[i][0] - np.median(po_col2[i][0])]
        phase_col_normalized2 += [po_col2[i][1] - np.median(po_col2[i][1])]

    Freq2 = np.array(freq_range2[0]) / 1e9

    X2, Y2 = np.meshgrid(Current2, Freq2)
    Z2 = np.array(np.transpose(logmag_col_normalized2))


# plt.pcolormesh(X, Y, Z,  vmin=-1., vmax=1., cmap='RdBu_r')
# plt.pcolormesh(X, Y, Z,  vmin=-5., vmax=5.)
plt.pcolormesh(X1, Y1, Z1, vmin=-10., vmax=10.,
               cmap='bwr', shading='auto')

plt.pcolormesh(X2, Y2, Z2, vmin=-2., vmax=2.,
               cmap='bwr', shading='auto')

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim([-1.4, -1.0])
plt.ylim([0, 7.5])
plt.xlabel("Current (mA)", fontsize=18)
plt.ylabel("Frequency (GHz)", fontsize=18)
# plt.colorbar()

# Click on the points on screen to define an approximation line for interpolation
def onclick(event):
    print('[%f, %f],' % (event.xdata, event.ydata))

# cid = fig.canvas.mpl_connect('button_press_event', onclick)


if 1:
    if 0:   # almost
        E_j = 5.401
        E_c = 1.228
        E_l = 1.483
        I_o = (0.8) * 1e-3  # One period
    if 0:
        E_l = 1.7235325179374792
        E_c = 1.1445558175924666
        E_j = 5.750693975688356
        I_o = (0.9) * 1e-3  # One period
    if 0:   # better than 0.9, 0.8
        E_l = 1.971333953257872
        E_c = 1.0752327771354064
        E_j = 6.112399561442605
        I_o = (1.0) * 1e-3  # One period
    if 1:   # better than 1.0
        E_l = 2.226396435599411
        E_c = 1.0161798344181885
        E_j = 6.483604035026594
        I_o = (1.1) * 1e-3  # One period
    if 0:   # better than 1.0
        E_l = 2.488336186528967
        E_c = 0.9650595627190067
        E_j = 6.86259657165787
        I_o = (1.2) * 1e-3  # One period
    current = np.linspace(-1.4, -1, 101) * 1e-3
    energy = np.zeros((len(current), 10))

    I_min = -1.109e-3  # half flux point current bias
    offset = (I_min - I_o / 2) / I_o  # 0 flux point

    phi_o = h / (2 * e)
    # flux1 = Current * phi_o / I_o
    # phi_ext1 = (flux1 / phi_o - offset) * 2 * np.pi

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
        # energy[idx, 2] = H.eigenenergies()[3] - H.eigenenergies()[0]
        # energy[idx, 3] = H.eigenenergies()[2] - H.eigenenergies()[1]
        # energy[idx, 1] = H.eigenenergies()[3] - H.eigenenergies()[2]
        # energy[idx, 2] = H.eigenenergies()[4] - H.eigenenergies()[3]
        # energy[idx, 3] = (H.eigenenergies()[3] - H.eigenenergies()[1])/2
        energy[idx, 4] = (H.eigenenergies()[2] - H.eigenenergies()[0])/2
        # energy[idx, 5] = (H.eigenenergies()[4] - H.eigenenergies()[2])/2
        # energy[idx, 1] = (H.eigenenergies()[3] - H.eigenenergies()[0])/3
        # energy[idx, 2] = (H.eigenenergies()[4] - H.eigenenergies()[1])/3

    cut = 400
    plt.plot(current * 1e3, energy[:, 0], 'k--', dashes=(3, 7), label='01')
    plt.plot(current * 1e3, energy[:, 1], 'y--', dashes=(3, 7), label='02')
    # plt.plot(current * 1e3, energy[:, 3], '--', dashes=(5, 5), label='12')
    plt.plot(current * 1e3, energy[:, 4], 'g--', dashes=(3, 7), label='02/2')
    # plt.plot(current * 1e3, energy[:, 5], '--', dashes=(5, 5), label='04/4')
    # plt.plot(current * 1e3, energy[:, 1], '--', dashes=(5, 5), label='23')
    # plt.plot(current * 1e3, energy[:, 2], '--', dashes=(5, 5), label='34')
    # plt.plot(current * 1e3, energy[:, 3], '--', dashes=(5, 5), label='13/2')
    # plt.plot(current * 1e3, energy[:, 5], '--', dashes=(5, 5), label='24/2')
    # plt.plot(current * 1e3, energy[:, 1], '--', dashes=(3, 7), label='03/3')
    # plt.plot(current * 1e3, energy[:, 2], '--', dashes=(3, 7), label='14/3')

    # plt.legend()

plt.title('Ej = 6.48 GHz, Ec = 1.02 GHz, El = 2.23 GHz', fontsize=18)
plt.show()
