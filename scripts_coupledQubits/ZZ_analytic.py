import numpy as np
from matplotlib import pyplot as plt


r = 6.0/6.41
na0 = -0.5743105321054214j
nb0 = 0.013663675268204414j
na1 = -0.5715754308834843j
nb1 = -0.03978926874384375j
ep = (na0 - r*na1)/(r*nb1 - nb0)
print (abs(ep))

freq2 = 4.484617233954527 #10-20
freq1 = 4.493589961059392 #11-21

# drive_amplitude_factor_array = np.linspace(0.2, 0.7, 51)*0.0524/0.5
# delta_omega_d_array = np.linspace(0.01,0.08,71)
# ZZ_rate = np.zeros((len(drive_amplitude_factor_array), len(delta_omega_d_array)))
#
# for a_idx, rabi0_1 in enumerate(drive_amplitude_factor_array):
#     for w_idx, detune1 in enumerate(delta_omega_d_array):
#         rabi0_2 = rabi0_1 * 6/6.41
#         detune2 = detune1 + (freq1 - freq2)
#         rabi_freq1 = np.sqrt(rabi0_1**2 + detune1**2)
#         AC_rate1 = 0.5*(rabi_freq1-detune1)
#         rabi_freq2 = np.sqrt(rabi0_2**2 + detune2**2)
#         AC_rate2 = 0.5*(rabi_freq2-detune2)
#         ZZ_rate[a_idx,w_idx] = abs(AC_rate2 - AC_rate1)
#
# X,Y = np.meshgrid(drive_amplitude_factor_array*500/0.0524, delta_omega_d_array + 4.49359)
# Z = ZZ_rate.transpose()*1e3 #MHz
# plt.pcolormesh(X,Y,Z, cmap = 'jet', vmax = 6.5,vmin=0.2)
# plt.ylabel('Drive freq (GHz)')
# plt.xlabel('Drive amplitude (mV)')
# plt.colorbar()
# plt.show()



