import numpy as np
from matplotlib import pyplot as plt
import scripts_singleQubit.plotting_settings
from scipy.optimize import curve_fit

fname = '/Users/longnguyen/Documents/tmp'
phase = np.load(fname+'_phaseZZ.npy')
print (phase.shape)

T_gate = 200
t_points = np.linspace(0, T_gate, 2*int(T_gate) + 1)

def linear_func(x,slope,offsetx, offsety):
    return slope*(x-offsetx)+offsety

ZZ_rate = np.zeros((51,61))
for a_idx in range(51):
    for w_idx in range(61):
        opt,cov = curve_fit(linear_func, ydata = np.unwrap(phase[a_idx,w_idx]), xdata = t_points)
        ZZ_rate[a_idx,w_idx]= abs(opt[0])
        # plt.plot(t_points, np.unwrap(phase[a_idx,w_idx]))
        # plt.plot(t_points, linear_func(t_points, *opt))
X,Y = np.meshgrid(np.linspace(0.025,0.075, 51)*500/0.05, np.linspace(0.03,0.06,61)+4.49)
Z = ZZ_rate.transpose()
plt.pcolormesh(X,Y,Z, cmap = 'jet_r')#, vmax = 0.0015)
# plt.colorbar()
plt.xlabel('Amplitude (mV)')
plt.ylabel('Frequency (GHz)')
plt.yticks([4.52,4.53, 4.54, 4.55])
plt.show()
