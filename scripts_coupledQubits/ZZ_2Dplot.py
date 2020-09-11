import numpy as np
from matplotlib import pyplot as plt
import scripts_singleQubit.plotting_settings
from scipy.optimize import curve_fit

width = 20
T_flat = 160
T_gate = T_flat+2*width
t_points = np.linspace(0, T_gate, 2*int(T_gate) + 1)

plt.figure(1, figsize = [9,6])
fname = '/Users/longnguyen/Documents/tmp'
phase = np.load(fname+'_phaseZZ_1121.npy')
print (phase.shape)
drive_amplitude_factor_array = np.linspace(0.2, 0.7, 51)*2*np.pi*0.0524/0.5
delta_omega_d_array = np.linspace(0.01,0.08,71)

def linear_func(x,slope,offsetx, offsety):
    return slope*(x-offsetx)+offsety

ZZ_rate = np.zeros((51,71))
for a_idx in range(51):
    for w_idx in range(71):
        opt,cov = curve_fit(linear_func, ydata = np.unwrap(phase[a_idx,w_idx]), xdata = t_points)
        ZZ_rate[a_idx,w_idx]= abs(opt[0])
        # plt.plot(t_points, np.unwrap(phase[a_idx,w_idx]))
        # plt.plot(t_points, linear_func(t_points, *opt))
X,Y = np.meshgrid(drive_amplitude_factor_array*500/0.0524/(2*np.pi), delta_omega_d_array + 4.4935)
Z = ZZ_rate.transpose()*1e3/(2*np.pi) #MHz
plt.pcolormesh(X,Y,Z+0.36, cmap = 'jet', vmax = 6.5,vmin=0.2)
plt.colorbar()
plt.xlabel('Amplitude (mV)')
plt.ylabel('Frequency (GHz)')
plt.yticks(np.linspace(4.51,4.56,6))
plt.ylim([4.51,4.56])

plt.figure(2, figsize = [10,7])
fname = '/Users/longnguyen/Documents/tmp'
phase = np.load(fname+'_phaseZZ_0131.npy')
print (phase.shape)

drive_amplitude_factor_array = np.linspace(0.2, 0.7, 51)*2*np.pi*0.091/0.5
delta_omega_d_array = np.linspace(0.03,0.085,56)

def linear_func(x,slope,offsetx, offsety):
    return slope*(x-offsetx)+offsety

ZZ_rate = np.zeros((51,56))
for a_idx in range(51):
    for w_idx in range(56):
        phase_fit = np.unwrap(phase[a_idx,w_idx])
        slope_guess = (phase_fit[-1] - phase_fit[0])/200
        guess = np.array([slope_guess,0,0])
        opt,cov = curve_fit(linear_func, ydata = phase_fit, xdata = t_points, p0 = guess)
        ZZ_rate[a_idx,w_idx]= abs(opt[0])
        # plt.plot(t_points, np.unwrap(phase[a_idx,w_idx]))
        # plt.plot(t_points, linear_func(t_points, *opt))
X,Y = np.meshgrid(drive_amplitude_factor_array*500/0.091/(2*np.pi), delta_omega_d_array + 6.6)
Z = ZZ_rate.transpose()*1e3/(2*np.pi) #MHz
plt.pcolormesh(X,Y,Z+0.36, cmap = 'jet', vmax = 13,vmin=1)
plt.colorbar()
plt.xlabel('Amplitude (mV)')
plt.ylabel('Frequency (GHz)')
plt.ylim([6.63,6.68])

plt.show()
