import os
import sys
script_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, script_folder)
from figure_properties import new_figure
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO


AR_inv = np.linspace(0, 1, 100)

data = """0	1
0.044093178	0.920187793
0.096505824	0.812206573
0.141430948	0.72143975
0.178868552	0.649452269
0.218801997	0.58372457
0.252911814	0.524256651
0.301996672	0.458528951
0.352745424	0.394366197
0.425124792	0.319248826
0.493344426	0.250391236
0.585690516	0.17370892
0.639767055	0.136150235
0.698003328	0.097026604
0.777870216	0.053208138
0.863560732	0.021909233
0.937603993	0.007824726
1	0"""

# Convert the string to a NumPy array
array = np.loadtxt(StringIO(data))

l_function = np.polyfit(array[:, 0], array[:, 1], 4)
lambda_c = np.polyval(l_function, AR_inv)


fig, ax = new_figure(latex_style=True)
scale = 0.6
fig.set_size_inches(6.4*scale, 4.8*scale)
ax.plot(AR_inv,lambda_c, lw=3, zorder =4, color='k')
ax.set_xlabel(r'$AR^{-1} = A_1/A_2$')
ax.set_ylabel(r'$\lambda = p_{loss}/(1/2 \rho V_1^2)$')
ax.grid(zorder=0)
fig.tight_layout(pad=0.5)

AR = np.linspace(1,4,100)
CP_I = 1 - (1/AR)**2
fig, ax = new_figure(latex_style=True)
scale = 0.75
fig.set_size_inches(6.4*scale, 4.8*scale)
lambda_c = np.polyval(l_function, 1/AR)
CP = CP_I - lambda_c
ax.plot(AR,CP_I, lw=3, zorder =4, color='b', label=r'Ideal $(C_{PR,I})$')
ax.plot(AR,CP, lw=3, zorder =4, color='r', label=r'Real $(C_{PR})$')
ax.set_xlabel(r'Area Ratio, $AR$')
ax.set_ylabel(r'$C_{PR}$')
ax.set_ylim([0,1])
ax.set_xticks(np.arange(1,4.1,0.5))
ax.set_yticks(np.arange(0, 1.1, 0.2))
ax.legend(draggable=True)
fig.tight_layout(pad=0.5)

nD = CP/CP_I
fig, ax = new_figure(latex_style=True)
scale = 0.6
fig.set_size_inches(6.4*scale, 4.8*scale)
ax.plot(AR,nD, lw=3, zorder =4, color='g')
ax.set_xlabel(r'Area Ratio, $AR$')
ax.set_ylabel(r'Diffuser Efficiency, $\eta_{D}$')
ax.set_xticks(np.arange(1,4.1,0.5))
ax.set_yticks(np.arange(0.4, 1.1, 0.1))
fig.tight_layout(pad=0.5)
plt.show()