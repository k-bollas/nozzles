import numpy as np
import matplotlib.pyplot as plt
from figure_properties import new_figure
from matplotlib.ticker import AutoMinorLocator
from labellines import *
import matplotlib as mpl

M1_range = np.arange(1.0, 3.01, 0.1)
fig, ax = new_figure(latex_style=True)
fig.set_size_inches(6.4, 4.8+2)
cmap = mpl.colormaps['Blues']

# Take colors at regular intervals spanning the colormap.
colors = cmap(np.linspace(0.4, 1, len(M1_range)))

for i,M1 in enumerate(M1_range):
    beta = np.linspace(0, 90, 100)

    theta = np.rad2deg(np.arctan(2/np.tan(np.deg2rad(beta))*(M1**2*np.sin(np.deg2rad(beta))**2-1)/(M1**2*(1.4+np.cos(np.deg2rad(2*beta)))+2)))


    line1, = ax.plot(theta,beta, lw=1.5, color=colors[i], zorder=3, label=f'{M1:1.1f}')
    ax.set_xlim([0, 35])
    ax.set_ylim([10,90])
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(10))  


ax.grid(zorder=0, which='major')
ax.grid(zorder=0, which='minor', lw=1, ls=':')
ax.set_xticks(np.arange(0, 35.1, 5))
ax.set_xticklabels([r'{:1.0f}$^\circ$'.format(x) for x in ax.get_xticks()])
ax.set_yticks(np.arange(0, 90.1, 10))
ax.set_yticklabels([r'{:1.0f}$^\circ$'.format(x) for x in ax.get_yticks()])
ax.set_xlabel(r'Angle $\theta$')
ax.set_ylabel(r'Angle $\beta$')
labelLines(ax.get_lines(), zorder=5.5)
norm = mpl.colors.Normalize(vmin=1.0, vmax=3.0)

fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax, orientation='horizontal', label='Mach Number')

# mplcursors.cursor(hover=True)

fig.tight_layout(pad=0.3)
# fig.savefig('theta_beta_equation.pgf')
plt.show()