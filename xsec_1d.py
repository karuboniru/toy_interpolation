import math
import numpy as np
import matplotlib.pyplot as plt

import sys
import os

sys.path.append('build')

import hkkm_interpolation

f = hkkm_interpolation.genie_spline('/var/home/yan/neutrino/spline/full/genie_spline_3_0_2/G00_00b_00_000/3.04.02-routine_validation_01-xsec_vA/total_xsec.root')

def func(E): 
    xsec = f.get_cross_section(14, 1000060120, E, "tot_cc")
    return xsec / E

x = np.linspace(0.1, 100, 3000)
y = np.array([func(logE) for logE in x])

plt.plot(x, y)
plt.xlabel('logE')
plt.ylabel('xsec / E')
plt.title('Plot of 1D Function')
plt.savefig('plot_xsec.eps')