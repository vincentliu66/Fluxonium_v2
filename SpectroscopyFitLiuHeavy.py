import numpy as np
from matplotlib import pyplot as plt
import sys
# sys.path.append('C:\Program Files (x86)\Labber\Script')
# import Labber
from qutip import *
from scipy.optimize import curve_fit

# import h5py

#####################################################################################
######################################Data###########################################
#####################################################################################

#############################################################################################

clicked_data1 = np.array([
    [-1.344436, 6.822484],
    [-1.326218, 6.689089],
    [-1.312807, 6.587709],
    [-1.304204, 6.502337],
    [-1.293071, 6.390285],
    [-1.279155, 6.235548],
    [-1.268275, 6.064802],
    [-1.259672, 5.931408],
    [-1.253852, 5.808685],
    [-1.244996, 5.621932],
    [-1.242719, 5.568575],
    [-1.227537, 5.168391],
    [-1.218934, 4.922945],
    [-1.173308, 3.092360],
    [-1.167680, 2.858407],
    [-1.164010, 2.718035],
    [-1.161319, 2.577662],
    [-1.158383, 2.443139],
    [-1.155447, 2.320314],
    [-1.153001, 2.203337],
    [-1.150799, 2.115604],
    [-1.148107, 2.004476],
    [-1.146150, 1.899197],
    [-1.142725, 1.770523],
    [-1.141012, 1.700336],
    [-1.139055, 1.606755],
    [-1.136119, 1.524871],
    [-1.133183, 1.396197],
    [-1.131226, 1.296767],
    [-1.127556, 1.185639],
    [-1.123886, 1.056964],
    [-1.121194, 0.951685],
    [-1.117524, 0.863952],
    [-1.097462, 0.916592],
    [-1.093303, 1.056964],
    [-1.090367, 1.173941],
    [-1.086697, 1.296767],
    [-1.084006, 1.384499],
    [-1.081559, 1.495627],
    [-1.079112, 1.618453],
    [-1.076666, 1.676941],
    [-1.075442, 1.752976],
    [-1.073485, 1.834860],
    [-1.071528, 1.899197],
    [-1.069081, 2.022023],
    [-1.067858, 2.086360],
    [-1.065411, 2.179941],
    [-1.063454, 2.255976],
    [-1.061986, 2.320314],
    [-1.059050, 2.466535],
    [-1.056359, 2.577662],
    [-1.054157, 2.665395],
    [-1.051465, 2.776523],
    [-1.049019, 2.911046],
    [-1.046327, 3.016325],
    [-1.043881, 3.109907],
    [-1.043147, 3.145000],

])

clicked_data2 = np.array([
    [0.472339, 8.243658],
    [0.486734, 7.891407],
    [0.505363, 7.603202],
    [0.542621, 7.250952],
    [0.568589, 6.984095],
    [0.583266, 6.770610],
    [0.684032, 6.823981],
    [0.717339, 7.186906],
    [0.759113, 7.581854],
    [0.794960, 8.297029]
])
current1 = clicked_data1[:, 0] * 1e-3  # In A
freq1 = clicked_data1[:, 1]  # in GHz

# current2 = clicked_data2[:,0]*1e-3 #In A
# freq2 = clicked_data2[:,1] #in GHz
current2 = []
freq2 = []

current = np.concatenate([current1, current2], axis=0)
freq = np.concatenate([freq1, freq2], axis=0)
# current = current1
# freq = freq1
plt.plot(current * 1e3, freq, 'o')  # plot mA
# plt.plot(current*1e3-1.023, freq, 'o') #plot mA
#####################################################################################
######################################Fit###########################################
#####################################################################################
# Define constants
e = 1.602e-19  # Fundamental charge
h = 6.62e-34  # Placnk's constant
phi_o = h / (2 * e)  # Flux quantum

N = 30
# E_l_guess = 1.629
# E_c_guess = 1.219
# E_j_guess = 7.6
# I_o = 1.023e-3
# offset = (0.633e-3 - I_o / 2) / I_o

E_l_guess = 2.5
E_c_guess = 0.96
E_j_guess = 6.8
I_o = 1.1 * 1e-3  # One period
I_min = -1.109e-3  # half flux point current bias
offset = (I_min - I_o / 2) / I_o  # 0 flux point

guess = ([E_l_guess, E_c_guess, E_j_guess])


def trans_energy(current, E_l, E_c, E_j):
    energy1 = np.zeros(len(current1))
    energy2 = np.zeros(len(current2))

    flux1 = current1 * phi_o / I_o
    phi_ext1 = (flux1 / phi_o - offset) * 2 * np.pi
    a = tensor(destroy(N))

    phi = (a + a.dag()) * (8.0 * E_c / E_l) ** (0.25) / np.sqrt(2.0)
    na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2.0)
    print('len(current1)=', len(current1))
    for idx in range(len(current1)):
        # print('idx=', idx)
        ope = 1.0j * (phi - phi_ext1[idx])
        # print('Here1')
        H = 4.0 * E_c * na ** 2 + 0.5 * E_l * phi ** 2 - 0.5 * E_j * (ope.expm() + (-ope).expm())
        # print('Here2')
        energy1[idx] = H.eigenenergies()[1] - H.eigenenergies()[0]

    # flux2 = current2 * phi_o / I_o
    # phi_ext2 = (flux2 / phi_o - offset) * 2 * np.pi
    # a = tensor(destroy(N))
    # phi = (a + a.dag()) * (8.0 * E_c / E_l) ** (0.25) / np.sqrt(2.0)
    # na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2.0)
    # for idx in range(len(current2)):
    #     ope = 1.0j * (phi - phi_ext2[idx])
    #     H = 4.0 * E_c * na ** 2.0 + 0.5 * E_l * phi ** 2.0 - 0.5 * E_j * (ope.expm() + (-ope).expm())
    #     energy2[idx] = H.eigenenergies()[2] - H.eigenenergies()[0]

    # return np.concatenate([energy1, energy2], axis=0)
    return energy1


opt, cov = curve_fit(trans_energy, current, freq, guess)

E_l_fit, E_c_fit, E_j_fit = opt
parameters_fit = {"E_l": E_l_fit, "E_c": E_c_fit, "E_j": E_j_fit}
for x, y in parameters_fit.items():
    print("{}={}".format(x, y))
# print ('E_l=' + str(E_l_fit) + ', E_c=' + str(E_c_fit) + ', E_j=' + str(E_j_fit) +
#        '\n' + 'A=' + str(A_fit) + ', offset='+ str(offset_fit))


############################################################################################################
# E_l, E_c, E_j = E_l_guess, E_c_guess, E_j_guess
E_l, E_c, E_j = E_l_fit, E_c_fit, E_j_fit
# current = np.linspace(-0.6, 1, 101) * 1e-3
current = np.linspace(-1.5, -0.75, 101) * 1e-3
energy = np.zeros((len(current), 10))

flux = current * phi_o / I_o
phi_ext = (flux / phi_o - offset) * 2 * np.pi
a = tensor(destroy(N))
phi = (a + a.dag()) * (8.0 * E_c / E_l) ** (0.25) / np.sqrt(2.0)
na = 1.0j * (a.dag() - a) * (E_l / (8 * E_c)) ** (0.25) / np.sqrt(2.0)
for idx in range(len(current)):
    ope = 1.0j * (phi - phi_ext[idx])
    H = 4.0 * E_c * na ** 2 + 0.5 * E_l * phi ** 2 - 0.5 * E_j * (ope.expm() + (-ope).expm())
    energy[idx, 0] = H.eigenenergies()[1] - H.eigenenergies()[0]
    energy[idx, 1] = H.eigenenergies()[2] - H.eigenenergies()[0]
    # energy[idx, 2] = H.eigenenergies()[3] - H.eigenenergies()[0]
    energy[idx, 2] = H.eigenenergies()[2] - H.eigenenergies()[1]

cut = 400
plt.plot(current * 1e3, energy[:, 0], '--')
# plt.plot(current * 1e3, energy[:, 1], '--')
# plt.plot(current*1e3, trans_energy(current,*guess),'s')
plt.show()
