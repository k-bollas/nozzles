import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

M1 = 1
M2 = 2.
A1 = 2
gamma = 1.4
A2 = (1/M2**2*(2/(gamma+1)*(1+(gamma-1)/2*M2**2))**((gamma+1)/(gamma-1)))*A1
mu1 = np.rad2deg(np.arcsin(1/M1))
mu2 = np.rad2deg(np.arcsin(1/M2))

nu_1 = np.rad2deg(np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt((gamma-1)/(gamma+1)*(M1**2-1))) - np.arctan(np.sqrt(M1**2-1)))
nu_2 = np.rad2deg(np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt((gamma-1)/(gamma+1)*(M2**2-1))) - np.arctan(np.sqrt(M2**2-1)))
theta_total = nu_2 - nu_1
theta_0 = theta_total/2
L = (A2/2 - A1/2)/np.tan(np.deg2rad(theta_0))

line_wall_plus = np.polyfit([0, L], [A1/2, A2/2], 1)

f, ax = plt.subplots()

line = np.zeros(2)
for mu in np.linspace(mu1, mu2, 20):
    line[0] = np.tan(np.deg2rad(mu))
    line[1] = -1 - line[0]*0
    y_zero = 0
    x_zero = (y_zero - line[1])/line[0]
    A = np.array([[-line_wall_plus[0], 1], [-line[0], 1]])
    b = np.array([[line_wall_plus[1]],[line[1]]])
    x_wall, y_wall = np.matmul(np.linalg.inv(A), b)
    x_wall, y_wall = float(x_wall), float(y_wall)
    if y_wall > A2/2:
        y_wall = A2/2
        x_wall = (y_wall - line[1])/line[0]
    # ax.plot([0, x_zero], [-1, y_zero], color='deepskyblue', linewidth = 0.5)
    # ax.plot([0, x_zero], [1, -y_zero], color='deepskyblue', linewidth = 0.5)
    # ax.plot([0, x_wall], [-1, y_wall], color='deepskyblue', linewidth = 0.5)
    # ax.plot([0, x_wall], [1, -y_wall], color='deepskyblue', linewidth = 0.5)

theta = theta_0
beta_fun = lambda x: np.tan(np.deg2rad(theta)) - 2/np.tan(np.deg2rad(x))*(M2**2*np.sin(np.deg2rad(x))**2 - 1)/(M2**2*(gamma + np.cos(np.deg2rad(2*x))) + 2)
beta_shock = fsolve(func= beta_fun, x0= 10, xtol=1e-5)[0]

shock_line = np.zeros(2)
shock_line[0] = np.tan(np.deg2rad(beta_shock - theta_0))
shock_line[1] = -A2/2 - shock_line[0]*L
y_zero_shock = 0
x_zero_shock = (y_zero_shock - shock_line[1])/shock_line[0]

ax.plot([L, x_zero_shock], [-A2/2, y_zero_shock], color='dodgerblue', linewidth = 3, linestyle='solid')
ax.plot([L, x_zero_shock], [A2/2, y_zero_shock], color='dodgerblue', linewidth = 3, linestyle='solid')

L_after = x_zero_shock + 1
ax.plot([0,0], [-1,1], color='black', linestyle = '--')
ax.plot([L,L],[-A2/2, A2/2], color='black', linestyle = '--')
ax.plot([0, L], [1, A2/2], color='blue', linewidth = 2.5)
ax.plot([0, L], [-1, -A2/2], color='blue', linewidth = 2.5)
ax.plot([-1, 0], [1, 1], color='blue', linewidth = 2.5)
ax.plot([L, L_after], [A2/2, A2/2], color='blue', linewidth = 2.5)
ax.plot([-1, 0], [-1, -1], color='blue', linewidth = 2.5)
ax.plot([L, L_after], [-A2/2, -A2/2], color='blue', linewidth = 2.5)
plt.show()

Mn2 = M2*np.sin(np.deg2rad(beta_shock))
Mn3 = np.sqrt(((1 + (gamma - 1)/2*Mn2**2))/(gamma*Mn2**2 - (gamma-1)/2))
M3 = Mn3/np.sin(np.deg2rad(beta_shock - theta))
Cp, R = 1008, 8314/29
# s_fun = lambda x: Cp*np.log((1 + 2*gamma/(gamma + 1)*(x**2 - 1))*(2 + (gamma - 1)*x**2)/((gamma + 1)*x**2))- R*np.log(1 + 2*gamma/(gamma + 1)*(x**2 - 1))
T_fun = lambda x: (1 + 2*gamma/(gamma+1)*(x**2 - 1))*(2 + (gamma-1)*x**2)/((gamma+1)*x**2)
P_fun = lambda x: 1 + 2*gamma/(gamma + 1)*(x**2-1)
s_fun = lambda x: Cp*np.log((1 + 2*gamma/(gamma+1)*(x**2 - 1))*(2 + (gamma-1)*x**2)/((gamma+1)*x**2)) - R*np.log(1 + 2*gamma/(gamma + 1)*(x**2-1))
Delta_s = s_fun(Mn2)

print(' Mach number after shock(s) : M_real_exit = {:.4f}'.format(M3))
print(' Entropy increase : {:.2f} J/kg-K'.format(Delta_s))

ax.set_xlim([-1, max(L_after, x_zero_shock)])
ax.set_ylim([-A2/2*1.02, A2/2*1.02])
ax.set_yticks([])
ax.set_xticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

line_wall_plus1 = np.zeros(2)
line_wall_plus1[0] = np.tan(np.deg2rad(theta_0))
line_wall_plus1[1] = A1/2 - line_wall_plus1[0]*0
f1, ax1 = plt.subplots()
for mu in np.linspace(mu1, mu2, 20):
    line[0] = np.tan(np.deg2rad(mu))
    line[1] = -1 - line[0]*0
    y_zero = 0
    x_zero = (y_zero - line[1])/line[0]
    A = np.array([[-line_wall_plus1[0], 1], [-line[0], 1]])
    b = np.array([[line_wall_plus1[1]],[line[1]]])
    x_wall, y_wall = np.matmul(np.linalg.inv(A), b)
    x_wall, y_wall = float(x_wall), float(y_wall)
    if y_wall > A2/2:
        y_wall = A2/2
        x_wall = (y_wall - line[1])/line[0]
    ax1.plot([0, x_zero], [-1, y_zero], color='deepskyblue', linewidth = 0.5)
    ax1.plot([0, x_zero], [1, -y_zero], color='deepskyblue', linewidth = 0.5)
    ax1.plot([0, x_wall], [-1, y_wall], color='deepskyblue', linewidth = 0.5)
    ax1.plot([0, x_wall], [1, -y_wall], color='deepskyblue', linewidth = 0.5)

L1 = x_wall
ax1.plot([0,0], [-1,1], color='black', linestyle = '--')
ax1.plot([0, L1], [1, y_wall], color='blue', linewidth = 2.5)
ax1.plot([0, L1], [-1, -y_wall], color='blue', linewidth = 2.5)

theta_1 = theta_0/2
line_wall_plus2 = np.zeros(2)
line_wall_plus2[0] = np.tan(np.deg2rad(theta_1))
line_wall_plus2[1] = y_wall - line_wall_plus2[0]*L1
L2 = (A2/2 - line_wall_plus2[1])/line_wall_plus2[0]
ax1.plot([L1, L2], [y_wall, A2/2], color='blue', linewidth = 2.5)
ax1.plot([L1, L2], [-y_wall, -A2/2], color='blue', linewidth = 2.5)
ax1.plot([L2,L2],[-A2/2, A2/2], color='black', linestyle = '--')
ax1.plot([-1, 0], [1, 1], color='blue', linewidth = 2.5)
ax1.plot([-1, 0], [-1, -1], color='blue', linewidth = 2.5)

theta = theta_1
beta_fun = lambda x: np.tan(np.deg2rad(theta)) - 2/np.tan(np.deg2rad(x))*(M2**2*np.sin(np.deg2rad(x))**2 - 1)/(M2**2*(gamma + np.cos(np.deg2rad(2*x))) + 2)
beta_shock1 = fsolve(func= beta_fun, x0= 10, xtol=1e-5)[0]
shock_line1 = np.zeros(2)
shock_line1[0] = np.tan(np.deg2rad(beta_shock1 - theta_1))
shock_line1[1] = -A2/2 - shock_line1[0]*L
y_zero_shock1 = 0
x_zero_shock1 = (y_zero_shock1 - shock_line1[1])/shock_line1[0]
ax1.plot([L1, x_zero_shock1], [-y_wall, -y_zero_shock1], color='dodgerblue', linewidth = 3, linestyle='solid')
ax1.plot([L1, x_zero_shock1], [y_wall, y_zero_shock1], color='dodgerblue', linewidth = 3, linestyle='solid')
Mn2 = M2*np.sin(np.deg2rad(beta_shock1))

Mn3 = np.sqrt(((1 + (gamma - 1)/2*Mn2**2))/(gamma*Mn2**2 - (gamma-1)/2))
M3 = Mn3/np.sin(np.deg2rad(beta_shock1 - theta_1))
Delta_s1 = s_fun(Mn2)
beta_fun = lambda x: np.tan(np.deg2rad(theta)) - 2/np.tan(np.deg2rad(x))*(M3**2*np.sin(np.deg2rad(x))**2 - 1)/(M3**2*(gamma + np.cos(np.deg2rad(2*x))) + 2)
beta_shock2 = fsolve(func= beta_fun, x0= 10, xtol=1e-5)[0]
shock_line2 = np.zeros(2)
shock_line2[0] = np.tan(np.deg2rad(beta_shock2 - theta_1))
shock_line2[1] = -A2/2 - shock_line2[0]*L2
y_zero_shock2 = 0
x_zero_shock2 = (y_zero_shock2 - shock_line2[1])/shock_line2[0]
ax1.plot([L2, x_zero_shock2], [-A2/2, -y_zero_shock2], color='dodgerblue', linewidth = 3, linestyle='solid')
ax1.plot([L2, x_zero_shock2], [A2/2, y_zero_shock2], color='dodgerblue', linewidth = 3, linestyle='solid')
L_after = x_zero_shock2 + 1
Mn3 = M3*np.sin(np.deg2rad(beta_shock2))
Mn4 = np.sqrt(((1 + (gamma - 1)/2*Mn3**2))/(gamma*Mn3**2 - (gamma-1)/2))
M4 = Mn4/np.sin(np.deg2rad(beta_shock2 - theta_1))
Delta_s2 = s_fun(Mn3)
Delta_s = Delta_s1 + Delta_s2


ax1.plot([L2, L_after], [A2/2, A2/2], color='blue', linewidth = 2.5)
ax1.plot([L2, L_after], [-A2/2, -A2/2], color='blue', linewidth = 2.5)
ax1.set_xlim([-1, max(L_after, x_zero_shock)])
ax1.set_ylim([-A2/2*1.02, A2/2*1.02])
ax1.set_yticks([])
ax1.set_xticks([])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)

print('--'*20)
print(' Mach number after first shock : M_real_exit (2) = {:.4f}'.format(M3))
print(' β1 = {:.4f} deg, β1 - θ1 = {:.4f} deg'.format(beta_shock1, beta_shock1-theta_1))
print(' Entropy increase : Δs1 {:.2f} J/kg-K'.format(Delta_s1))
print(' β2 = {:.4f} deg, β2 - θ1 = {:.4f} deg'.format(beta_shock2, beta_shock2-theta_1))
print(' Mach number after first shock : M_real_exit = {:.4f}'.format(M4))
print(' Entropy increase : Δs2 {:.2f} J/kg-K'.format(Delta_s2))
print(' Total Entropy increase : Δs {:.2f} J/kg-K'.format(Delta_s))

plt.show()