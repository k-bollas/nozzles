import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import matplotlib.ticker as ticker
from figure_properties import new_figure

A_A_max = 10.0
M1, = fsolve(lambda M: np.sqrt(1/M**2*(2/2.4*(1+0.2*M**2))**(2.4/0.4)) - A_A_max, x0 = 0.01, xtol=1e-5 )
M2, = fsolve(lambda M: np.sqrt(1/M**2*(2/2.4*(1+0.2*M**2))**(2.4/0.4)) - A_A_max, x0 = 5, xtol=1e-5 )

print(M1, M2)
M = np.logspace(np.log10(M1), np.log10(M2),200)
A_A = np.sqrt(1/M**2*(2/2.4*(1+0.2*M**2))**(2.4/0.4))

M_scatter_1 = np.arange(1,3.51, 0.5)
A_A_scatter = np.sqrt(1/M_scatter_1**2*(2/2.4*(1+0.2*M_scatter_1**2))**(2.4/0.4))
M_scatter_2 = np.array([fsolve(lambda M: np.sqrt(1/M**2*(2/2.4*(1+0.2*M**2))**(2.4/0.4)) - y, x0=0.05, xtol=1e-5) for y in A_A_scatter])
M_scatter = np.append(M_scatter_1, M_scatter_2)
A_A_scatter = np.append(A_A_scatter, A_A_scatter)

fig, ax = new_figure(latex_style=True)
scale = 1.2
fig.set_size_inches(6.4*scale, 4.8*scale+0.1)
ax.plot(A_A, M, color='b', lw=3.5, zorder=0)
ax.scatter(A_A_scatter, M_scatter, marker='o', color='b')
ax.set_yscale('log')
M_labels = [0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0]
ax.set_yticks(M_labels)
ax.set_yticklabels([f"{x:1.2f}" if np.round(x,1) != x else f"{x:1.1f}" for x in M_labels])
ax.set_xticks(np.arange(1,10.1))
ax.grid(which='both', lw=1.25, color=(0.4,0.4,0.4,1), zorder=1)
ax.set_xlabel(r'Area Ratio, $A^{*}/A$')
ax.set_ylabel(r'Mach Number, $M$')
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))    # Minor ticks every 1 on y-axis

ax.axhline(1.0, color='k', lw=2.5, ls='-', zorder=5)
fig.tight_layout(pad=0.5)
# fig.savefig('area_mach_distr.pgf')

fig, ax = new_figure(latex_style=True)
scale = 0.6
fig.set_size_inches(6.4*scale, 4.8*scale+0.1)
ax.plot(A_A, M, color='b', lw=3.5, zorder=0)
ax.scatter(A_A_scatter, M_scatter, marker='o', color='b')
ax.set_yscale('log')
M_labels = [0.1, 0.5, 1.0, 4.0]
ax.set_yticks(M_labels)
ax.set_yticklabels([f"{x:1.2f}" if np.round(x,1) != x else f"{x:1.1f}" for x in M_labels])
ax.set_xticks(np.arange(1,10.1))
ax.grid(which='both', lw=0.5, color=(0.4,0.4,0.4,1), zorder=1)
ax.set_xlabel(r'Area Ratio, $A^{*}/A$')
ax.set_ylabel(r'Mach Number, $M$')
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))    # Minor ticks every 1 on y-axis

ax.axhline(1.0, color='k', lw=2.5, ls='-', zorder=5)
fig.tight_layout(pad=0.5)
# fig.savefig('area_mach_distr.pgf')

plt.show()