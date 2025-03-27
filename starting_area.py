import os
import sys
script_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, script_folder)
from figure_properties import new_figure
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

poly = [0.2435, - 0.113, + 0.8701]
M = np.linspace(1,3,200)
As_Ao = np.polyval(poly, M)

fig, ax = new_figure(latex_style=True)
scale = 0.6
fig.set_size_inches(6.4*scale, 4.8*scale)
ax.plot(M, As_Ao, '-b', lw=3.5, zorder=4)
ax.set_xlabel(r'$M_{FL}$ (-)')
ax.set_ylabel(r'$\frac{A_{s}}{A_{o}}$', rotation=0, fontsize=25, labelpad=20)
ax.set_xticks([1,2,3])
ax.set_yticks([1,2,3])
fig.tight_layout(pad=0.5)
plt.show()