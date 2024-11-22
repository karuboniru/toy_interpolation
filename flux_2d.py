import hkkm_interpolation
import math
import numpy as np
import matplotlib.pyplot as plt

import sys
import os

sys.path.append('build')

from mpl_toolkits.mplot3d import Axes3D

interpolator = hkkm_interpolation.hkkm_2d('/var/home/yan/code/toy_interpolation/data/honda-2d.txt')
def func(logE, costh):
    E = 10**logE
    flux = interpolator.get_flux(E, costh, 14)
    return math.log10(flux)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.linspace(-1, 4, 300)
y = np.linspace(-1, 1, 300)

X, Y = np.meshgrid(x, y)

Z = np.array([[func(logE, costh) for logE in x] for costh in y])

# Plot the surface.
ax.plot_surface(X, Y, Z, cmap='viridis')

ax.set_xlabel('logE')
ax.set_ylabel('cos(theta)')
ax.set_zlabel('log10(flux)')
plt.title('Plot of 2D Function')

# Save the figure to plot.eps
plt.savefig('plot.eps')