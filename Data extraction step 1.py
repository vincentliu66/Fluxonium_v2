# Plot intensity data

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
rc('text', usetex=False)

#Enter directory and name of measurement
directory = 'G:\Projects\Fluxonium & qubits\Data\\2016_01\\13'
measurement = 'S21_Phase2tones_ZNB_0&n40_BW100Hz_SMB_n20dBm_6to8GHz_YOKO_n7ton10V'
path_data = directory + '\\' + measurement + '_Phase.csv'
path_freq = directory + '\\' + measurement + '_Freq.csv'
path_vol = directory + '\\' + measurement + '_Voltage.csv'

RawData = np.genfromtxt(path_data, delimiter =',')
Freq = np.genfromtxt(path_freq, delimiter =',')/1e9
V = np.genfromtxt(path_vol, delimiter =',')

Z = RawData.transpose()
          
#Optional: calculate differential
Z_diff = np.diff(Z.transpose()) 
Z_diff = Z_diff.transpose()

#Plot the intensity map
X, Y = np.meshgrid(V,Freq)
fig=plt.figure()
plt.pcolormesh(X, Y, Z, cmap=cm.BuGn, vmin =-1 , vmax =1)
plt.title(measurement)
plt.xlabel("Voltage (V)")
plt.ylabel("Frequency (GHz)")
plt.colorbar()

#Click on the points on screen to define an approximation line for interpolation
def onclick(event):
    print '[%f, %f],'%(event.xdata, event.ydata)
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()