import os
import sys
script_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, script_folder)
from figure_properties import new_figure
import matplotlib.pyplot as plt
import numpy as np

AR = np.linspace(1,5,200)
CPR = 1 - 1 / AR**2
AR_p = np.arange(1,5.1,0.5)
CPR_p = 1 - 1/AR_p**2

fig, ax = new_figure(latex_style=True)
scale = 0.7
fig.set_size_inches(6.4*scale, 4.8*scale)

ax.plot(AR, CPR, lw=4, color='b', zorder=4)
ax.scatter(AR_p, CPR_p, s=40, color='b', zorder=5)
ax.grid()
ax.set_xlabel(r'Diffuser area ratio, $AR$')
ax.set_ylabel(r'Ideal $C_{PR}$')
plt.tight_layout(pad=0.5)
plt.show()