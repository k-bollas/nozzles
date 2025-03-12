from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sc
import matplotlib.patches as patc
from figure_properties import new_figure
from NozzleMOC__ import NozMOC
import pandas as pd

M1 = 0.2
M_son = 1
M2 = 2
g = 1.4

p2 = 101325
p1_p_son = (1 + (g-1)/2*M_son**2)**(g/(g-1))/(1 + (g-1)/2*M1**2)**(g/(g-1))
p_son_p2 = (1 + (g-1)/2*M2**2)**(g/(g-1))/(1 + (g-1)/2*M_son**2)**(g/(g-1))
p_son = p_son_p2*p2
p1 = p_son*p1_p_son
p0 = (1 + (g-1)/2*M1**2)**(g/(g-1))*p1
p0_p_crit = (2/(g+1))**(g/(g-1))

A_son = 1
A1_A_son = np.sqrt(1/M1**2*(2/(g+1)*(1+(g-1)/2*M1**2))**((g+1)/(g-1)))
A2_A_son = np.sqrt(1/M2**2*(2/(g+1)*(1+(g-1)/2*M2**2))**((g+1)/(g-1)))
A1 = A_son*A1_A_son
A2 = A_son*A2_A_son

print('Exit Pressure : p2 = {:.4f} bar'.format(p2/1e5))
print('Sonic Pressure : p* = {:.4f} bar'.format(p_son/1e5))
print('Inlet Pressure : p1 = {:.4f} bar'.format(p1/1e5))
print('Total Pressure : p0 = {:.4f} bar'.format(p0/1e5))
print('')
print('Sonic Area : A* = {:.2f} m2'.format(A_son))
print('Inlet Area : A1 = {:.2f} m2'.format(A1))
print('Exit Area : A2 = {:.2f} m2'.format(A2))
R = A1 - A_son
L = R
L_list_sub = np.linspace(0,L,40)
A_list_sub, M_list_sub, p_p0_list_sub = [], [], []
conv_wall = np.polyfit([0,L],[A1/2,A_son/2], 1)

for Li in L_list_sub:
    theta = np.arccos((Li-R)/R)
    A_list_sub.append(A1 - R*np.sin(theta))
    M_list_sub.append(sc.fsolve(func= lambda x: 1/x**2*(2/(g+1)*(1+(g-1)/2*x**2))**((g+1)/(g-1)) - (A_list_sub[-1]/A_son)**2, x0= np.random.rand(1)*0.1, xtol=1e-5)[0])
    p0_p = (1 + (g-1)/2*M_list_sub[-1]**2)**(g/(g-1))
    p_p0_list_sub.append(1/p0_p)

div_wall = np.polyfit([L,2*L],[A_son/2,A2/2], 1)
L_list_sup = np.linspace(L,2*L,40)
A_list_sup, M_list_sup, M_list_sub_lim, p_p0_list_sup, p_p0_list_sup_lim = [], [], [], [], []
for Li in L_list_sup:
    A_list_sup.append(np.polyval(div_wall, Li)*2)
    M_list_sup.append(sc.fsolve(func= lambda x: 1/x**2*(2/(g+1)*(1+(g-1)/2*x**2))**((g+1)/(g-1)) - (A_list_sup[-1]/A_son)**2, x0= 1.5, xtol=1e-5)[0])
    M_list_sub_lim.append(sc.fsolve(func= lambda x: 1/x**2*(2/(g+1)*(1+(g-1)/2*x**2))**((g+1)/(g-1)) - (A_list_sup[-1]/A_son)**2, x0= np.random.rand(1)*0.1, xtol=1e-5)[0])
    p0_p = (1 + (g-1)/2*M_list_sup[-1]**2)**(g/(g-1))
    p_p0_list_sup.append(1/p0_p)
    p0_p = (1 + (g-1)/2*M_list_sub_lim[-1]**2)**(g/(g-1))
    p_p0_list_sup_lim.append(1/p0_p)

f1, ax1 = new_figure(latex_style=True)
scale = 0.8
f1.set_size_inches(6.4*scale, 4.8*scale)
plt1 = ax1.plot(L_list_sub, p_p0_list_sub, lw=4, color='black', label='Pressure distribution', zorder=4)
ax1.plot(L_list_sup, p_p0_list_sup, lw=4, color='black', zorder=4)
ax1.plot([0,L],[p0_p_crit, p0_p_crit], color='black', linestyle='--', linewidth=1, zorder=3)
ax1.plot([L,L],[0, p0_p_crit], color='black', linestyle='--', linewidth=1, zorder=3)
ax1.set_yticks([round(x,y) for x,y in zip([p_p0_list_sup[-1], p0_p_crit, 1],[3,3,0])])
ax1.set_yticklabels([str(round(x,y)) for x,y in zip([0, p0_p_crit, 1],[0,3,0])])
ax1.set_ylabel(r'Relative pressure : $p/p0$')
ax1.set_xlim([0,2*L*1.05])
ax1.set_ylim([min(p_p0_list_sup)*0.95, 1])
ax1.set_xticks([0,L,2*L])
ax1.set_xticklabels(['Inlet', 'Throat', 'Exit'])
ax1.grid(zorder=0)
# f1.savefig()

ax2 = ax1.twinx()
plt2 = ax2.plot(L_list_sub, M_list_sub, color='gray', lw=4, linestyle='--', label='Mach Number distribution', zorder=4)
ax2.plot(L_list_sup, M_list_sup, color='gray', lw=4, linestyle='--', zorder=4)
ax2.plot([L,2*L],[1, 1], color='black', linestyle='--', linewidth=1, zorder=3)
ax2.plot([L,L],[0, 1], color='black', linestyle='--', linewidth=1, zorder=3)
ax2.set_ylabel('Mach Number')
ax2.set_xlim([0,2*L*1.00])
ax2.set_ylim([0, max(M_list_sup)+0.18])
ax2.set_yticks([0, 1, max(M_list_sup)])
lns = plt1 + plt2
labs = [l.get_label() for l in lns]
legend = ax2.legend(lns, labs, loc="lower center", bbox_to_anchor=(0.5, 1.03), facecolor='white')
# ax1.set_ylim([ax1.get_ylim()[0], ax2.get_ylim()[1]/max(M_list_sup)])
ax2.grid(zorder=0)
f1.tight_layout(pad=0.3)
f1.savefig('pressure_mach_distr.eps', format='eps', dpi=600)


f2, (ax1, ax3, ax4) = plt.subplots(1,3)
f2.set_figheight(5)
f2.set_figwidth(15)
plt1 = ax1.plot(L_list_sub, p_p0_list_sub, color='black', label='p/p0 ratio')
ax1.plot(L_list_sup, p_p0_list_sup, color='blue')
ax1.plot([0,L],[p0_p_crit, p0_p_crit], color='black', linestyle='--')
ax1.plot([L,L],[0, p0_p_crit], color='black', linestyle='--')
ax1.set_yticks([round(x,y) for x,y in zip([p_p0_list_sup[-1], p0_p_crit, 1],[3,3,0])])
ax1.set_ylabel('p/p0')
ax1.set_xlim([0,2*L*1.05])
ax1.set_ylim([min(p_p0_list_sup)*0.95, 1])
ax1.set_xticks([0,L,2*L])
ax1.set_xticklabels(['Inlet', 'Throat', 'Exit'])

ax1.legend(['Nominal case'])

ax2 = ax1.twinx()
plt2 = ax2.plot(L_list_sub, M_list_sub, color='green', label='Mach Number')
ax2.plot(L_list_sup, M_list_sup, color='green')
ax2.plot([L,2*L],[1, 1], color='black', linestyle='--')
ax2.plot([L,L],[0, 1], color='black', linestyle='--')
ax2.set_ylabel('Mach Number')
ax2.set_xlim([0,2*L*1.00])
ax2.set_ylim([0, max(M_list_sup)*1.275])
ax2.set_yticks([0, 1, 1.8])
lns = plt1 + plt2
labs = [l.get_label() for l in lns]
legend = ax1.legend(lns, labs, loc="upper right", bbox_to_anchor=(1.03, 1.03), edgecolor="black")
legend.get_frame().set_alpha(0)
legend.get_frame().set_facecolor((1, 1, 1, 1))


ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.get_xaxis().set_ticks([])
ax3.get_yaxis().set_ticks([])

ax4.plot(L_list_sub, p_p0_list_sub, color='black')
ax4.plot(L_list_sup, p_p0_list_sup, color='black')
ax4.plot(L_list_sup, p_p0_list_sup_lim, color='black', linestyle='--')
ax4.plot([0,L],[p0_p_crit, p0_p_crit], color='black', linestyle='--')
ax4.plot([L,L],[0, p0_p_crit], color='black', linestyle='--')
ax4.set_yticks([round(x,3) for x in [1, p0_p_crit, 0]])
ax4.set_yticklabels(["1", r'$\left(\frac{2}{\gamma + 1}\right)**{\frac{\gamma}{\gamma -1}}$', "0"])
ax4.set_ylabel(r'p/$p_{t}$')
ax4.set_xlim([0,2*L*1.05])
ax4.set_ylim([0, 1])
ax4.set_xticks([0,L,2*L])
ax4.set_xticklabels(['Inlet', 'Throat', 'Exit'])

f2.suptitle('Nominal Case (p_exit = p_atm)', fontsize=16)
ax3.legend(['Nominal case'])
# plt.savefig('Figure2.png', transparent=True)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Case p_exit > p_atm
p2 = 0.96*p0
M2 = np.sqrt(2/(g-1)*((p0/p2)**((g-1)/g) -1))
A2_A_son = np.sqrt(1/M2**2*(2/(g+1)*(1+(g-1)/2*M2**2))**((g+1)/(g-1)))
A_son = A2/A2_A_son
p_p0_ = 1/((1 + (g-1)/2*M2**2)**(g/(g-1)))

M_list_sub2, p_p0_list_sub2 = [], []
for i in range(0,len(A_list_sub)):
    M_list_sub2.append(sc.fsolve(func= lambda x: 1/x**2*(2/(g+1)*(1+(g-1)/2*x**2))**((g+1)/(g-1)) - (A_list_sub[i]/A_son)**2, x0= 0.1, xtol=1e-5)[0])
    p0_p = (1 + (g-1)/2*M_list_sub2[-1]**2)**(g/(g-1))
    p_p0_list_sub2.append(1/p0_p)

M_list_sup2, p_p0_list_sup2 = [], []
for i in range(0,len(A_list_sup)):
    M_list_sup2.append(sc.fsolve(func= lambda x: 1/x**2*(2/(g+1)*(1+(g-1)/2*x**2))**((g+1)/(g-1)) - (A_list_sup[i]/A_son)**2, x0= 0.1, xtol=1e-5)[0])
    p0_p = (1 + (g-1)/2*M_list_sup2[i]**2)**(g/(g-1))
    p_p0_list_sup2.append(1/p0_p)

print('------------------------------------------------------------------------------------------------------------------------------')
print('New Sonic Area : A* new = {:.4f} -'.format(A_son))
print('Exit Mach Number : M1 = {:.4f} -'.format(M2))
p_p0_list_sub1, p_p0_list_sup1 = p_p0_list_sub, p_p0_list_sup

f3, (ax1, ax2, ax3) = plt.subplots(1,3)
f3.set_figheight(5)
f3.set_figwidth(15)
print(pd.DataFrame(np.transpose(np.array([L_list_sub, p_p0_list_sub2, M_list_sup2]))))
# ax1.plot(L_list_sub, p_p0_list_sub, color='blue', linestyle = '--')
# ax1.plot(L_list_sup, p_p0_list_sup, color='blue', linestyle = '--')
ax1.plot(L_list_sub, p_p0_list_sub2, color='blue')
ax1.plot(L_list_sup, p_p0_list_sup2, color='blue')

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.get_xaxis().set_ticks([])
ax2.get_yaxis().set_ticks([])

ax1.plot([0,L],[p0_p_crit, p0_p_crit], color='black', linestyle='--')
ax1.plot([L,L],[0, p0_p_crit], color='black', linestyle='--')
ax1.set_xticklabels([])
ax1.set_yticks([round(x,3) for x in [p2/p0]])
# ax1.set_yticks([round(x,3) for x in [p_p0_list_sup[-1], p0_p_crit, p2/p0, 1]])
# ax1.set_ylabel('p/p0')
ax1.grid(axis='y')
# ax1.set_xticks([0,L,2*L])
# ax1.set_xticklabels(['Inlet', 'Throat', 'Exit'])
ax1.set_xlim([0,2*L*1.05])
ax1.set_ylim([min(p_p0_list_sup)*0.95, 1])

ax3.plot(L_list_sub, M_list_sub, color='green', linestyle = '--')
ax3.plot(L_list_sup, M_list_sup, color='green', linestyle = '--')
ax3.plot(L_list_sub, M_list_sub, color='green')
ax3.plot(L_list_sup, M_list_sub_lim, color='green')
ax3.set_xticklabels([])
ax3.set_yticks([round(x,3) for x in [M_list_sub_lim[-1]]])
# ax3.set_yticks([round(x,3) for x in [0, 1, M2, 1.8]])
ax3.grid(axis='y')
# ax3.set_ylabel('M')
ax3.set_xlim([0,2*L*1.05])
ax3.set_ylim([0, max(M_list_sup)*1.05])
# ax3.set_xticks([0,L,2*L])
# ax3.set_xticklabels(['Inlet', 'Throat', 'Exit'])
# plt.savefig('Figure3.png', transparent=True)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Case p_exit > p_crit (M* -> 1 -)
M2 = 2
p2 = 101325
M_exit = np.sqrt((1 + (g-1)/2*M2**2)/(g*M2**2 - (g-1)/2))
p_exit = (1 + 2*g/(g+1)*(M2**2 - 1))*p2
p_exit_p0 = p_exit/p0


f4, ax4 = new_figure()

ax4.plot(L_list_sub, M_list_sub, color='black', linewidth = 3)
ax4.plot(L_list_sup, M_list_sup, color='black', linewidth = 3)
ax4.plot([L_list_sup[-1], L_list_sup[-1]], [M_list_sup[-1], 1.234287053], color='black', linewidth = 3)
ax4.set_xticks([0,L,2*L])
ax4.set_xticklabels(['Inlet', 'Throat', 'Exit'])
ax4.set_yticks([0,1.234287053, 1, 2])
ax4.set_yticklabels(['0', r'$M_{exit}$', '1', r'$M_{des,exit}$'])
ax4.set_title('Mach Number distribution in C-D Nozzle', fontsize=18)
ax4.set_ylabel('Mach number', fontsize=18)
ax4.grid(linewidth=3)
ax4.set_xlim([0,2*L*1.05])
ax4.set_ylim([0, max(M_list_sup)*1.05])
ax4.tick_params(axis='both', which='major', labelsize=18)
plt.tight_layout()
plt.savefig('Mach_shock_oblique.pgf')


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Case p_exit > p_crit (M* -> 1 -)
M_before_shock = 1.56
M_after_shock = np.sqrt((1 + (g-1)/2*M_before_shock**2)/(g*M_before_shock**2 - (g-1)/2))
A_son = 1
A_beforeshock_A_son = np.sqrt(1/M_before_shock**2*(2/(g+1)*(1+(g-1)/2*M_before_shock**2))**((g+1)/(g-1)))
A_shock = A_beforeshock_A_son*A_son
L_shock = sc.fsolve(func= lambda x: np.polyval(div_wall, x)*2 - A_beforeshock_A_son, x0 = L, xtol=1e5)[0]
p_ratio = 1 + 2*g/(g+1)*(M_before_shock**2 - 1)
T_ratio = p_ratio*(2 + (g-1)*M_before_shock**2)/((g+1)*M_before_shock**2)
Cp, R = 1008, 8314/29
Ds = Cp*np.log(T_ratio) - R*np.log(p_ratio)
p0_ratio = np.exp(-Ds/R)

A_aftershock_A_son = np.sqrt(1/M_after_shock**2*(2/(g+1)*(1+(g-1)/2*M_after_shock**2))**((g+1)/(g-1)))
A_son_aftershock = A_shock/A_aftershock_A_son
A_list_sup3_beforeshock = np.linspace(A_list_sup[0], A_shock, 20)
L_list_sup3, M_list_sup3, p_p0_list_sup3 = [], [], []
for i in range(0, len(A_list_sup3_beforeshock)):
    L_list_sup3.append(sc.fsolve(func= lambda x: np.polyval(div_wall, x)*2 - A_list_sup3_beforeshock[i], x0 = L, xtol=1e5)[0])
    M_list_sup3.append(sc.fsolve(func= lambda x: 1/x**2*(2/(g+1)*(1+(g-1)/2*x**2))**((g+1)/(g-1)) - (A_list_sup3_beforeshock[i]/A_son)**2, x0= 1.5, xtol=1e-5)[0])
    p0_p = (1 + (g-1)/2*M_list_sup3[-1]**2)**(g/(g-1))
    p_p0_list_sup3.append(1/p0_p)

A_list_sup3_aftershock = np.linspace(A_shock, A_list_sup[-1], 20)
for i in range(0, len(A_list_sup3_beforeshock)):
    L_list_sup3.append(sc.fsolve(func= lambda x: np.polyval(div_wall, x)*2 - A_list_sup3_aftershock[i], x0 = L, xtol=1e5)[0])
    M_list_sup3.append(sc.fsolve(func= lambda x: 1/x**2*(2/(g+1)*(1+(g-1)/2*x**2))**((g+1)/(g-1)) - (A_list_sup3_aftershock[i]/A_son_aftershock)**2, x0= np.random.rand(1)*0.1, xtol=1e-5)[0])
    p0_p = (1 + (g-1)/2*M_list_sup3[-1]**2)**(g/(g-1))*(1/p0_ratio)
    p_p0_list_sup3.append(1/p0_p)

f5, ax4 = new_figure()

ax4.plot(L_list_sub, M_list_sub, color='black', linestyle = '--')
ax4.plot(L_list_sup, M_list_sup, color='black', linestyle = '--')
ax4.plot(L_list_sub, M_list_sub, color='black', linewidth = 3)
ax4.plot(L_list_sup3, M_list_sup3, color='black', linewidth = 3)
ax4.set_xticks([0,L,2*L])
ax4.set_xticklabels(['Inlet', 'Throat', 'Exit'])
ax4.set_yticks([0,M_list_sup3[-1], 1, 2])
ax4.set_yticklabels(['0', r'$M_{exit}$', '1', r'$M_{exit,des}$'])
ax4.set_title('Mach Number distribution in C-D Nozzle', fontsize=18)
ax4.set_ylabel('Mach number', fontsize=18)
ax4.grid(linewidth=3)
ax4.set_xlim([0,2*L*1.05])
ax4.set_ylim([0, max(M_list_sup)*1.05])
ax4.tick_params(axis='both', which='major', labelsize=18)
plt.tight_layout()
# plt.savefig('Mach_shock.pgf')


f5, ax5 = new_figure(latex_style=True, box=True)
ax5.plot(L_list_sub, p_p0_list_sub, linewidth=2, color='black')
ax5.plot(L_list_sup, p_p0_list_sup, linewidth=2, color='black')
ax5.plot(L_list_sup, p_p0_list_sup_lim, color='black', linewidth=2, linestyle=(0, (5, 1)) )
ax5.plot(L_list_sub, p_p0_list_sub2, linewidth=1.5, color='black', linestyle=(5, (10, 3)))
ax5.plot(L_list_sup, p_p0_list_sup2, linewidth=1.5, color='black', linestyle=(5, (10, 3)))
ax5.plot(L_list_sub, p_p0_list_sub,  linewidth=1.5, color='black', linestyle=(0, (1, 1)))
ax5.plot(L_list_sup3, p_p0_list_sup3,  linewidth=1.5, color='black', linestyle=(0, (1, 1)))
ax5.plot([L_list_sup[-1], L_list_sup[-1]], [p_p0_list_sup[-1], 0.45], linewidth=1.5, color='black', linestyle=(0, (3, 1, 1, 1)))
ax5.text(L_list_sup[-8], p_p0_list_sup2[-8], "A", ha="center", va="center", bbox=dict(boxstyle="round", fc=(1., 1, 1),))
ax5.text(L_list_sup[15], p_p0_list_sup_lim[15], "B", ha="center", va="center", bbox=dict(boxstyle="round", fc=(1., 1, 1),))
ax5.text(1.03*L_list_sup[-1], p_p0_list_sup3[-1], "C", ha="center", va="center", bbox=dict(boxstyle="round", fc=(1., 1, 1),))
ax5.plot([L_list_sup[-1], 1.1*L_list_sup[-1]], [0.45, 0.55], linewidth=1, color='black')
ax5.text(1.1*L_list_sup[-1], 0.55, "D", ha="center", va="center", bbox=dict(boxstyle="round", fc=(1., 1, 1),))
ax5.plot([L_list_sup[-1], 1.1*L_list_sup[-1]], [0.3, 0.3], linewidth=1, color='black')
ax5.text(1.1*L_list_sup[-1], 0.3, "E", ha="center", va="center", bbox=dict(boxstyle="round", fc=(1., 1, 1),))
furter_expansion = np.polyfit([L_list_sup[-1], 1.05*L_list_sup[-1], 1.2*L_list_sup[-1]], [p_p0_list_sup[-1], 0.8*p_p0_list_sup[-1], 0.6*p_p0_list_sup[-1]], 2)
ax5.plot(np.linspace(L_list_sup[-1], 1.2*L_list_sup[-1], 50), np.polyval(furter_expansion, np.linspace(L_list_sup[-1], 1.2*L_list_sup[-1], 50)), linewidth=1.5, color='black', linestyle='--')
ax5.plot([L_list_sup[-1], 1.1*L_list_sup[-1]], [p_p0_list_sup[-1], 0.2], linewidth=1, color='black')
ax5.text(1.1*L_list_sup[-1], 0.2, "F", ha="center", va="center", bbox=dict(boxstyle="round", fc=(1., 1, 1),))
ax5.plot([1.05*L_list_sup[-1], 0.9*L_list_sup[-1]], [0.8*p_p0_list_sup[-1], 0.6*p_p0_list_sup[-1]], linewidth=1, color='black')
ax5.text(0.9*L_list_sup[-1], 0.6*p_p0_list_sup[-1], "G", ha="center", va="center", bbox=dict(boxstyle="round", fc=(1., 1, 1),))

ax5.grid(linewidth= 3)
ax5.set_yticks([round(x,3) for x in [1, p0_p_crit, 0]])
ax5.set_yticklabels(["1", r'$\left(\frac{2}{\gamma + 1}\right)^{\frac{\gamma}{\gamma -1}}$', "0"])
ax5.set_ylabel(r'p/$p_{t}$', fontsize=14)
ax5.set_xlim([0,1.2*L_list_sup[-1]])
ax5.set_ylim([0, 1])
ax5.set_xticks([0,L,2*L])
ax5.set_xticklabels(['Inlet', 'Throat', 'Exit'])

ax5.set_title('Pressure distribution in C-D Nozzle', fontsize=14)
ax5_ = ax5.twinx()
ax5_.set_yticks([p_p0_list_sup_lim[-1], 0.45, p_p0_list_sup[-1]])
ax5_.set_yticklabels([r'$p_{crit,1}$', r'$p_{crit,2}$', r'$p_{exit,des}$'])
ax5_.grid(linewidth=1.5)
ax5.tick_params(axis='both', which='major', labelsize=14)
ax5_.tick_params(axis='both', which='major', labelsize=14)

plt.tight_layout()
plt.savefig('pressure_nozzle.pgf')

# ======================================================================================================================================================================================
# 
# NOZZLE GEOMETRY
# 
# ======================================================================================================================================================================================

L = 8
a_mat = np.array([[L**3, L**2], [3*L**2, 2*L]])
c_mat = np.array([[A_son - A1], [0]])
b_mat = np.matmul(np.linalg.inv(a_mat), c_mat)

x_sub = np.linspace(0, L, 100)
y_sub = np.array([b_mat[0] *x**3 + b_mat[1]*x**2 + A1 for x in x_sub])

x0, y0 = L, 1
x_super, y_super, _ = NozMOC(2.4, 30, g, x0, y0, y0)
f6, ax6 = new_figure(box=True)
scale = 0.8
f6.set_size_inches(6.4*scale, 4.8*scale)
ax6.plot(x_sub, y_sub, lw=4, color='blue')
ax6.plot(x_sub, -y_sub, lw=4, color='blue')
ax6.plot(x_super, y_super, lw=4, color='blue')
ax6.plot(x_super, -y_super, lw=4, color='blue')
ax6.plot([x0, x0], [-y0, y0], color='blue', linewidth=1.5, linestyle='--')
ax6.set_xlabel('Nozzle Length (X)')
ax6.set_ylabel('Nozzle Length (Y)')
ax6.set_xticks([*np.linspace(0, L, 4), *np.linspace(L, x_super[-1], 4)])
ax6.set_xticklabels(['Inlet', *['' for x in range(3)], 'Throat', *['' for x in range(2)], 'Exit'])
ax6.set_yticklabels('')
# ax6.grid()
f6.tight_layout(pad=0.3)
# f6.savefig('cd_nozzle.eps', dpi=600, format='eps')
f6.savefig('cd_nozzle.svg', dpi=600, format='svg')


f7, ax7 = new_figure()
ax7.plot(x_sub, y_sub, linewidth = 4, color='blue')
ax7.plot(x_sub, -y_sub, linewidth = 4, color='blue')
ax7.plot(x_super, y_super, linewidth = 4, color='blue')
ax7.plot(x_super, -y_super, linewidth = 4, color='blue')
ax7.plot([x0, x0], [-y0, y0], color='blue', linewidth=1.5, linestyle='--')
ax7.set_xlabel('Nozzle Length (X)', fontsize=18)
ax7.set_ylabel('Nozzle Length (Y)', fontsize=18)
ax7.set_xticks([*np.linspace(0, L, 4), *np.linspace(L, x_super[-1], 4)])
ax7.set_xticklabels(['Inlet', *['' for x in range(3)], 'Throat', *['' for x in range(2)], 'Exit'])
ax7.set_yticklabels('')
ax7.set_title('Convergent-Divergent Nozzle', fontsize=18)
ax7.tick_params(axis='both', which='major', labelsize=18)
ax7.plot([x_super[int(len(x_super)/2)], x_super[int(len(x_super)/2)]], [-y_super[int(len(x_super)/2)], y_super[int(len(x_super)/2)]], linewidth=1.5, color='black')
ax7.plot(np.array([x_super[int(len(x_super)/2)], x_super[int(len(x_super)/2)]])*1.01, [-y_super[int(len(x_super)/2)], y_super[int(len(x_super)/2)]], linewidth=1.5, color='black')
ax7.text(x_super[int(len(x_super)/2)] - 0.5, 0, r'$M > 1$', ha="right", va="center", size=20)
ax7.text(x_super[int(len(x_super)/2)] + 0.5, 0, r'$M < 1$', ha="left", va="center", size=20)
ax7.quiver(x_super[int(len(x_super)/2)] - 3.5, -0.55 , 3.5, 0, color='black', scale=21)
ax7.quiver(x_super[int(len(x_super)/2)] + 1, -0.55 , 1.5, 0, color='black', scale=21)
ax7.plot([x_super[int(len(x_super)/2)], x_super[int(len(x_super)/2)] + 0.6], [y_super[int(len(x_super)/2)], y_super[int(len(x_super)/2)] + 0.5], linewidth=1, color='black')
ax7.text(x_super[int(len(x_super)/2)] + 0.5, y_super[int(len(x_super)/2)] + 0.6 , 'Shock Wave', ha="left", va="center", size=14)
plt.tight_layout()
# plt.savefig('CDnozzle_shockwave.pgf')

f8, ax8 = new_figure()
ax8.plot(x_sub, y_sub, linewidth = 4, color='blue')
ax8.plot(x_sub, -y_sub, linewidth = 4, color='blue')
ax8.plot(x_super, y_super, linewidth = 4, color='blue')
ax8.plot(x_super, -y_super, linewidth = 4, color='blue')
ax8.plot([x0, x0], [-y0, y0], color='blue', linewidth=1.5, linestyle='--')
ax8.set_xlabel('Nozzle Length (X)', fontsize=18)
ax8.set_ylabel('Nozzle Length (Y)', fontsize=18)
ax8.set_xticks([*np.linspace(0, L, 4), *np.linspace(L, x_super[-1], 4)])
ax8.set_xticklabels(['Inlet', *['' for x in range(3)], 'Throat', *['' for x in range(2)], 'Exit'])
ax8.set_yticklabels('')
ax8.set_title('Convergent-Divergent Nozzle', fontsize=18)
ax8.tick_params(axis='both', which='major', labelsize=18)
ax8.plot([x_super[-1], x_super[-1]], [-y_super[-1], y_super[-1]], linewidth=1.5, color='black')
ax8.plot(np.array([x_super[-1], x_super[-1]])*0.99, [-y_super[-1], y_super[-1]], linewidth=1.5, color='black')
ax8.text(x_super[-1] - 0.5, 0, r'$M > 1$', ha="right", va="center", size=20)
ax8.text(x_super[-1] + 0.5, 0, r'$M < 1$', ha="left", va="center", size=20)
ax8.quiver(x_super[-1] - 3.5, -0.55 , 3.5, 0, color='black', scale=21)
ax8.plot([x_super[-1], x_super[-1] + 0.6], [y_super[-1], y_super[-1] + 0.5], linewidth=1, color='black')
ax8.text(x_super[-1] - 0.5, y_super[-1] + 0.6 , 'Shock Wave', ha="left", va="center", size=14)
plt.tight_layout()
# plt.savefig('CDnozzle_shockwave_exit.pgf')

obl_shock = np.zeros(2)
obl_shock[0] = np.tan(np.deg2rad(40))
obl_shock[1] = (-y_super[-1]) - obl_shock[0]*(x_super[-1])

f9, ax9 = new_figure()
ax9.plot(x_sub, y_sub, linewidth = 4, color='blue')
ax9.plot(x_sub, -y_sub, linewidth = 4, color='blue')
ax9.plot(x_super, y_super, linewidth = 4, color='blue')
ax9.plot(x_super, -y_super, linewidth = 4, color='blue')
ax9.plot([x0, x0], [-y0, y0], color='blue', linewidth=1.5, linestyle='--')
ax9.set_xlabel('Nozzle Length (X)', fontsize=18)
ax9.set_ylabel('Nozzle Length (Y)', fontsize=18)
ax9.set_xticks([*np.linspace(0, L, 4), *np.linspace(L, x_super[-1], 4)])
ax9.set_xticklabels(['Inlet', *['' for x in range(3)], 'Throat', *['' for x in range(2)], 'Exit'])
ax9.set_yticklabels('')
ax9.set_title('Convergent-Divergent Nozzle', fontsize=18)
ax9.tick_params(axis='both', which='major', labelsize=18)
ax9.plot([x_super[-1], 1.2*x_super[-1]], np.polyval(obl_shock, [x_super[-1], 1.2*x_super[-1]]), linewidth=1.5, color='black')
obl_shock[1] = (-y_super[-1]) - obl_shock[0]*(x_super[-1]) + 0.15
ax9.plot(np.array([x_super[-1]-0.1, 1.2*x_super[-1]])*1, np.polyval(obl_shock, np.array([x_super[-1]-0.1, 1.2*x_super[-1]])*1), linewidth=1.5, color='black')
obl_shock[0] = -np.tan(np.deg2rad(40))
obl_shock[1] = (y_super[-1]) - obl_shock[0]*(x_super[-1])
ax9.plot([x_super[-1], 1.2*x_super[-1]], np.polyval(obl_shock, [x_super[-1], 1.2*x_super[-1]]), linewidth=1.5, color='black')
obl_shock[1] = (y_super[-1]) - obl_shock[0]*(x_super[-1]) - 0.15
ax9.plot(np.array([x_super[-1]-0.1, 1.2*x_super[-1]])*1, np.polyval(obl_shock, np.array([x_super[-1]-0.1, 1.2*x_super[-1]])*1), linewidth=1.5, color='black')
ax9.plot([x_super[-1], 1.2*x_super[-1]], [-y_super[-1], -y_super[-1]], linewidth=1.5, color='black', linestyle='--')
ax9.plot([x_super[-1], 1.2*x_super[-1]], [y_super[-1], y_super[-1]], linewidth=1.5, color='black', linestyle='--')
a1 = patc.Arc((x_super[-1], -y_super[-1]), width=0.1*x_super[-1], height=0.1*x_super[-1], angle=0.0, theta1=0.0, theta2=40)
ax9.add_patch(a1)
ax9.text(1.1*x_super[-1], -0.9*y_super[-1], r'$\beta$', ha="center", va="center", size=18)
ax9.text(x_super[-1], 0, r'$M > 1$', ha="right", va="center", size=20)
ax9.quiver(x_super[-1] - 3.5, -0.55 , 3.5, 0, color='black', scale=21)
obl_shock[1] = (y_super[-1]) - obl_shock[0]*(x_super[-1])
ax9.plot([1.1*x_super[-1], x_super[-1] + 0.6], [np.polyval(obl_shock, 1.1*x_super[-1]), y_super[-1] + 0.5], linewidth=1, color='black', linestyle=':')
ax9.text(x_super[-1] - 0.5, y_super[-1] + 0.6 , 'Oblique Shock Wave', ha="center", va="center", size=14)
plt.tight_layout()
# plt.savefig('CDnozzle_shockwave_oblique.pgf')

f10, ax10 = new_figure(latex_style=True)
ax10.plot(x_sub, y_sub, linewidth = 4, color='blue')
ax10.plot(x_sub, -y_sub, linewidth = 4, color='blue')
ax10.plot(x_super, y_super, linewidth = 4, color='blue')
ax10.plot(x_super, -y_super, linewidth = 4, color='blue')
ax10.plot([x0, x0], [-y0, y0], color='blue', linewidth=1.5, linestyle='--')
ax10.set_xlabel('Nozzle Length (X)', fontsize=18)
ax10.set_ylabel('Nozzle Length (Y)', fontsize=18)
ax10.set_xticks([*np.linspace(0, L, 4), *np.linspace(L, x_super[-1], 4)])
ax10.set_xticklabels(['Inlet', *['' for x in range(3)], 'Throat', *['' for x in range(2)], 'Exit'])
ax10.set_yticklabels('')
ax10.set_title('Convergent-Divergent Nozzle', fontsize=18)
ax10.tick_params(axis='both', which='major', labelsize=18)

angle_range = np.linspace(17,22,5)
for i,angle in enumerate(angle_range):
    exp_wave = np.zeros(2)
    exp_wave[0] = np.tan(np.deg2rad(angle))
    exp_wave[1] = (-y_super[-1]) - exp_wave[0]*(x_super[-1])
    extr_factor = 1.3
    ax10.plot([x_super[-1], extr_factor*x_super[-1]], np.polyval(exp_wave, [x_super[-1], extr_factor*x_super[-1]]), linewidth=0.5, color='black')
    exp_wave[0] = -np.tan(np.deg2rad(angle))
    exp_wave[1] = (y_super[-1]) - exp_wave[0]*(x_super[-1])
    ax10.plot([x_super[-1], extr_factor*x_super[-1]], np.polyval(exp_wave, [x_super[-1], extr_factor*x_super[-1]]), linewidth=0.5, color='black')
ax10.text(14.130566355312295 , -1.4462180181908204, r'$M_{1} > 1$', ha="left", va="center", size=15)
ax10.quiver(14.130566355312295 , -1.9009668648170157 , 2.5, 0, color='black', scale=21)
ax10.quiver(18.09182890252353 , -2.0146540764735645 , 2.5, -0.7, color='black', scale=21)
ax10.quiver(14.130566355312295 , 1.9009668648170157 , 2.5, 0, color='black', scale=21)
ax10.quiver(18.09182890252353 , 2.0146540764735645 , 2.5, 0.7, color='black', scale=21)
ax10.plot([18.967695289726855, 17.17425459212005], [1.5475452220982984, 2.7981045503203354], linewidth=1, color='black', linestyle=':')
ax10.text(17.17425459212005 , 3.025478973633433 , 'Expansion Waves', ha="center", va="center", size=14)
ax10.text(20.883238157660447 , -1.882018996207591 , r'$M_{2} > M_{1}$', ha="center", va="center", size=15, rotation=-13)

line_after_nozzle = np.zeros(2)
line_after_nozzle[0] = np.tan(np.deg2rad(-7.8))
line_after_nozzle[1] = (-y_super[-1]) - line_after_nozzle[0]*(x_super[-1])
ax10.plot([x_super[-1], 1.2*x_super[-1]], [-y_super[-1], np.polyval(line_after_nozzle, 1.2*x_super[-1])], linewidth=1.5, color='black', linestyle='--')
line_after_nozzle[0] = np.tan(np.deg2rad(7.8))
line_after_nozzle[1] = (y_super[-1]) - line_after_nozzle[0]*(x_super[-1])
ax10.plot([x_super[-1], 1.2*x_super[-1]], [y_super[-1], np.polyval(line_after_nozzle, 1.2*x_super[-1])], linewidth=1.5, color='black', linestyle='--')
plt.tight_layout()
f10.savefig('CDnozzle_expansion_waves.pgf')


def mouse_event(event):
    print('{} , {}'.format(event.xdata, event.ydata))

cid = f10.canvas.mpl_connect('button_press_event', mouse_event)


plt.show()
