import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

M_list = np.linspace(1.2,3,50)
theta = 26.43/2
gamma = 1.4
beta_shock = [fsolve(func= lambda x: np.tan(np.deg2rad(theta)) - 2/np.tan(np.deg2rad(x))*(y**2*np.sin(np.deg2rad(x))**2 - 1)/(y**2*(gamma + np.cos(np.deg2rad(2*x))) + 2), x0= theta, xtol=1e-5)[0] for y in M_list]
Mn2 = [x*np.sin(np.deg2rad(y)) for x,y in zip(M_list, beta_shock)]
Mn3 = [np.sqrt(((1 + (gamma - 1)/2*x**2))/(gamma*x**2 - (gamma-1)/2)) for x in Mn2]
M3 = [x/np.sin(np.deg2rad(y - theta)) for x,y in zip(Mn3, beta_shock)]


plt.figure()
plt.plot(M_list, beta_shock)
plt.figure()
plt.plot(M_list, M3)
plt.show()