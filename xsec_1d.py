import math
import numpy as np
import matplotlib.pyplot as plt

import sys
import os

sys.path.append('build')

import hkkm_interpolation

f = hkkm_interpolation.genie_spline('/var/home/yan/neutrino/spline/full/genie_spline_3_0_2/G18_10a_02_11b/3.04.02-routine_validation_01-xsec_vA/total_xsec.root')

def func(E): 
    return f.get_cross_section(12, 1000060120, E, "tot_cc") / E  + f.get_cross_section(12, 2212, E, "tot_cc") / E

x = np.logspace(-1, 2, 3000)
y = func(x)

plt.plot(np.log10(x), y)
plt.xlabel('logE')
plt.ylabel('xsec / E')
plt.title('Plot of 1D Function')
plt.show()