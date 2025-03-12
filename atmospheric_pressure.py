from ambiance.ambiance import Atmosphere
import numpy as np
import matplotlib.pyplot as plt
from figure_properties import new_figure

height = np.linspace(0, 50000, 100)
P_amb = np.array([Atmosphere(x).pressure for x in height])

scale = 0.6
fig, ax = new_figure(latex_style=True)
fig.set_size_inches(6.4*scale, 4.8*scale)
ax.plot(height/1e3, P_amb/1e3, lw=4, color='C2')
ax.set_xlabel(r'Altitude $\times 10^3$ (m)')
ax.set_ylabel('Ambient Pressure (kPa)')
ax.set_xticks(np.arange(0,51,10))
ax.set_yticks(np.arange(0,101,20))
plt.tight_layout(pad=0.3)
# plt.savefig('ambient_pressure_altitude.pgf')
plt.show()
