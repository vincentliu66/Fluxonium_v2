"""
CH Liu calculate the thermal occupation of fluxonium
Let's only consider the lowest-four energy levels
"""

import numpy as np

from scipy.constants import *

E0 = 0
E1 = h*5e9
E2 = h*(5+4.7)*1e9
E3 = h*(5+4.7+4.4)*1e9
E = np.array([E0, E1, E2, E3])
aef