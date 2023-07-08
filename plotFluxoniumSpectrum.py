"""
Chuan-Hong Liu, plot the spectrum of fluxonium
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import matplotlib.cm as cm

# datadir_1 = 'F:\Projects\2023-Fluxonium_Debug\PlotLightHeavyFluxoniums\Q3' \
#             '-current_sweep_pulsed_full_with_res_fitting_vs_flux_0.csv'

datadir_1 = 'Q3-current_sweep_pulsed_full_with_res_fitting_vs_flux_0.csv'
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

# print('freq_range=', np.array(freq_range[0]))


###

logmag_col = []
phase_col = []

logmag_col_normalized = []
phase_col_normalized = []

for i in range(np.shape(po_col)[0]):
    logmag_col += [po_col[i][0]]
    phase_col += [po_col[i][1]]

    logmag_col_normalized += [po_col[i][0] - np.median(po_col[i][0])]
    phase_col_normalized += [po_col[i][1] - np.median(po_col[i][1])]

# data.index[-1][0] = data.index[-1][0] - 4.95 * 1e-3

# fig, ax = plt.subplots(1, 2, figsize=(20, 7))

# np.flip( axis =1)

# Original
# im_amp = ax[0].imshow(np.transpose(logmag_col_normalized), cmap='viridis', interpolation=None, origin='lower',
#                       vmin=-0.2, vmax=0.3, aspect='auto',
#                       extent=(data.index[-1][0] * 1e3, data.index[0][0] * 1e3, freq_range[0][0], freq_range[0][-1]))
# im_phase = ax[1].imshow(np.flip(np.transpose(phase_col_normalized)), cmap='viridis', interpolation=None, origin='lower',
#                         vmin=-0.007, vmax=0.013, aspect='auto',
#                         extent=(data.index[-1][0] * 1e3, data.index[0][0] * 1e3, freq_range[0][0], freq_range[0][-1]))

# im_amp = ax[0].imshow(np.transpose(logmag_col_normalized), cmap='viridis', interpolation=None, origin='lower',
#                       vmin=-0.2, vmax=0.3, aspect='auto',
#                       extent=(data.index[0][0] * 1e3, data.index[-1][0] * 1e3, freq_range[0][0], freq_range[0][-1]))
# im_phase = ax[1].imshow(np.transpose(phase_col_normalized), cmap='viridis', interpolation=None, origin='lower',
#                         vmin=-0.007, vmax=0.013, aspect='auto',
#                         extent=(data.index[0][0] * 1e3, data.index[-1][0] * 1e3, freq_range[0][0], freq_range[0][-1]))

# ax[0].axhline(2.1*1e9, color = 'black')
# ax[0].axvline(4.95, color = 'black')
# fig.colorbar(im_amp, ax=ax[0])
# plt.show()


