import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import numpy as np
import scipy.special as sp

mpl.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

ax.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
circle1 = plt.Circle((0, 0), 1.5, color='lightblue')
ax.add_patch(circle1)

ax.set_xlim(0,5)
ax.set_ylim(0,3)
ax.set_xlabel('x')
ax.set_aspect('equal')

ax.annotate('(4), (5), (6), (7), (12)', xy=(math.sqrt(1.5**2/2), math.sqrt(1.5**2/2)),  xycoords='data',
            xytext=(3, 1.5), textcoords='data',
            arrowprops=dict(facecolor='black', shrink=0.05, width = 0.5, headwidth = 5),
            horizontalalignment='right', verticalalignment='top',
            )

plt.show()