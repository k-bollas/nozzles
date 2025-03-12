import numpy as np
import matplotlib.pyplot as plt
import os
import sys
script_folder = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_folder)
from figure_properties import new_figure
from matplotlib.ticker import AutoMinorLocator
from labellines import *
from matplotlib.ticker import AutoMinorLocator


M = np.arange(1,5,0.05)
nu = np.rad2deg(np.sqrt(2.4/0.4)*np.arctan(np.sqrt(0.4/2.4*(M**2-1))) - np.arctan(np.sqrt(M**2-1)))

fig, ax = new_figure(latex_style=True)
scale = 0.8
fig.set_size_inches(6.4*scale, 4.8*scale)
ax.plot(M,nu,'-b', lw=4, zorder=4)
ax.set_xlabel(r'Mach number, $M$')
ax.set_xticks(np.arange(1,5.1,0.5))
ax.set_ylabel(r'Prandlt-Meyer function $\nu(M)$ (deg)', y=0.45)
ax.set_yticks(np.arange(0,81,10))

ax.grid(zorder=0)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))  
ax.grid(which='minor', lw=1)
fig.tight_layout(pad=0.8)
fig.savefig('anlge_prandlt_meyer.pgf')
plt.show()