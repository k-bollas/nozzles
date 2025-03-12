import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sc
from figure_properties import new_figure

M1 = 0.05
M_son = 1
M2 = 1.8
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

fig1, ax1 = new_figure(latex_style=True)
scale = 0.8
fig1.set_size_inches(6.4*scale, 4.8*scale+0.5)
plt1 = ax1.plot(L_list_sub, p_p0_list_sub, lw=2.5, color='blue', label=r'$p/p_{t}$', zorder=4)
ax1.plot(L_list_sup, p_p0_list_sup, lw=2.5, color='blue', zorder=4)
ax1.plot([0,L],[p0_p_crit, p0_p_crit], color='blue', linestyle='--', linewidth=1, zorder=3)
ax1.plot([L,L],[0, p0_p_crit], color='blue', linestyle='--', linewidth=1, zorder=3)
ax1.set_yticks([round(x,y) for x,y in zip([p2/p1, p0_p_crit, 1],[3,3,0])])
ax1.set_yticklabels([str(round(x,y)) for x,y in zip([p2/p1, p0_p_crit, 1],[3,3,0])])
ax1.set_ylabel(r'$p/p_{t}$')
ax1.set_xlim([0,2*L*1.05])
ax1.set_ylim([p2/p1*0.25, 1])
ax1.set_xticks([0,L,2*L])
ax1.set_xticklabels(['Inlet', 'Throat', 'Exit'])
ax1.grid(zorder=0)

ax2 = ax1.twinx()
plt2 = ax2.plot(L_list_sub, M_list_sub, color='green', lw=2.5, linestyle='-', label='Mach Number', zorder=4)
ax2.plot(L_list_sup, M_list_sup, color='green', lw=2.5, linestyle='-', zorder=4)
ax2.plot([L,2*L],[1, 1], color='green', linestyle='--', linewidth=1, zorder=3)
ax2.plot([L,L],[0, 1], color='green', linestyle='--', linewidth=1, zorder=3)
ax2.set_ylabel('Mach Number')
ax2.set_xlim([0,2*L*1.00])
ax2.set_ylim([0, max(M_list_sup)+0.18])
ax2.set_yticks([0, 1, max(M_list_sup)])
lns = plt1 + plt2
labs = [l.get_label() for l in lns]
legend = ax2.legend(lns, labs, facecolor='white', loc='upper center', draggable=True)
# ax1.set_ylim([ax1.get_ylim()[0], ax2.get_ylim()[1]/max(M_list_sup)])
ax2.grid(zorder=0)
fig1.tight_layout(pad=0.5)

fig21, ax21 = new_figure(latex_style=True, box=True)
scale = 0.8
fig21.set_size_inches(6.4*scale, 4.8*scale+0.5)
ax21.plot(L_list_sub, p_p0_list_sub, lw=2.5, color='blue', linestyle='--', label=r'$p/p_{t}$', zorder=4)
ax21.plot(L_list_sup, p_p0_list_sup, lw=2.5, color='blue', linestyle='--', zorder=4)
ax21.plot([0,L],[p0_p_crit, p0_p_crit], color='blue', linestyle='--', linewidth=1, zorder=3)
ax21.plot([L,L],[0, p0_p_crit], color='blue', linestyle='--', linewidth=1, zorder=3)
ax21.set_yticks([round(x,y) for x,y in zip([p2/p1, p0_p_crit, 1],[3,3,0])])
ax21.set_yticklabels([str(round(x,y)) for x,y in zip([p2/p1, p0_p_crit, 1],[3,3,0])])
ax21.set_ylabel(r'$p/p_{t}$')
ax21.set_xticks([0,L,2*L])
ax21.set_xticklabels(['Inlet', 'Throat', 'Exit'])
ax21.set_xlim([0,2*L*1.05])
ax21.set_ylim([p2/p1*0.25, 1])
ax21.grid(zorder=0)
fig21.tight_layout(pad=0.5)

fig22, ax22 = new_figure(latex_style=True, box=True)
scale = 0.8
fig22.set_size_inches(6.4*scale, 4.8*scale+0.5)
plt2 = ax22.plot(L_list_sub, M_list_sub, color='green', lw=2.5, linestyle='--', label='Mach Number', zorder=4)
ax22.plot(L_list_sup, M_list_sup, color='green', lw=2.5, linestyle='--', zorder=4)
ax22.plot([L,2*L],[1, 1], color='green', linestyle='--', linewidth=1, zorder=3)
ax22.plot([L,L],[0, 1], color='green', linestyle='--', linewidth=1, zorder=3)
ax22.set_ylabel('Mach Number')
ax22.set_xlim([0,2*L*1.00])
ax22.set_ylim([0, max(M_list_sup)+0.18])
ax22.set_yticks([0, 1, max(M_list_sup)])
ax22.grid(zorder=0)
ax22.set_xticks([0,L,2*L])
ax22.set_xticklabels(['Inlet', 'Throat', 'Exit'])
fig22.tight_layout(pad=0.5)


# Case pe/p0 = 0.95
pp0 = 0.868
p2 = pp0*p0
M2 = np.sqrt(2/(g-1)*((p0/p2)**((g-1)/g) -1))
A2_A_son = np.sqrt(1/M2**2*(2/(g+1)*(1+(g-1)/2*M2**2))**((g+1)/(g-1)))
A_son = A2/A2_A_son
p_p0_ = 1/((1 + (g-1)/2*M2**2)**(g/(g-1)))

M_list_sub2, p_p0_list_sub2 = [], []
for i in range(0,len(A_list_sub)):
    M_list_sub2.append(sc.fsolve(func= lambda x: 1/x**2*(2/(g+1)*(1+(g-1)/2*x**2))**((g+1)/(g-1)) - (A_list_sub[i]/A_son)**2, x0=0.005, xtol=1e-5)[0])
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

fig31, ax31 = new_figure(latex_style=True, font_size=18)
scale = 0.8
fig31.set_size_inches(6.4*scale, 4.8*scale+0.5)
ax31.plot(L_list_sub, p_p0_list_sub2, color='blue', lw=2.5)
ax31.plot(L_list_sup, p_p0_list_sup2, color='blue', lw=2.5)
# ax31.plot(L_list_sub, p_p0_list_sub, lw=2.5, color='r', linestyle='--', label=r'$p/p_{t}$', zorder=4)
# ax31.plot(L_list_sup, p_p0_list_sup, lw=2.5, color='r', linestyle='--', zorder=4)

ax31.plot([0,2*L],[pp0, pp0], color='black', linestyle='--')
ax31.set_ylim(ax21.get_ylim())
ax31.set_xticks([0,L,2*L])
ax31.set_xticklabels([])
ax31.set_yticks([pp0])
ax31.set_xlim(ax21.get_xlim())
# ax31.set_xticklabels(['Inlet', 'Throat', 'Exit'], color='white')
fig31.tight_layout(pad=0.5)
fig31.savefig('fig31.svg', transparent=True)

fig32, ax32 = new_figure(latex_style=True, font_size=18)
scale = 0.8
fig32.set_size_inches(6.4*scale, 4.8*scale+0.5)

ax32.plot(L_list_sub, M_list_sub2, color='green', linestyle='-', lw=2.5)
ax32.plot(L_list_sup, M_list_sup2, color='green', linestyle='-', lw=2.5)
ax32.set_xlim(ax22.get_xlim())
ax32.set_ylim(ax22.get_ylim())
ax32.set_yticks([M_list_sup2[-1]])
ax32.grid(zorder=0)
ax32.set_xticks([0,L,2*L])
ax32.set_xticklabels([])
# ax31.set_xticklabels(['Inlet', 'Throat', 'Exit'], color='white')
fig32.tight_layout(pad=0.5)
fig32.savefig('fig32.svg', transparent=True)

# Normal shock

# Case p_exit > p_crit (M* -> 1 -)
M_before_shock = 1.8
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

fig41, ax41 = new_figure(latex_style=True, font_size=18)
scale = 0.8
fig41.set_size_inches(6.4*scale, 4.8*scale+0.5)
ax41.plot(L_list_sub, p_p0_list_sub2, color='blue', lw=2.5)
ax41.plot(L_list_sup3, p_p0_list_sup3, color='blue', lw=2.5)
# ax41.plot(L_list_sub, p_p0_list_sub, lw=2.5, color='r', linestyle='--', label=r'$p/p_{t}$', zorder=4)
# ax41.plot(L_list_sup, p_p0_list_sup, lw=2.5, color='r', linestyle='--', zorder=4)
ax41.plot([0,2*L],[p_p0_list_sup3[-1], p_p0_list_sup3[-1]], color='black', linestyle='--')
ax41.set_ylim(ax21.get_ylim())
ax41.set_xticks([0,L,2*L])
ax41.set_xticklabels([])
ax41.set_yticks([p_p0_list_sup3[-1]])
ax41.set_yticklabels([np.round(p_p0_list_sup3[-1],3)])
ax41.set_xlim(ax21.get_xlim())
# ax31.set_xticklabels(['Inlet', 'Throat', 'Exit'], color='white')
fig41.tight_layout(pad=0.5)
fig41.savefig('fig41.svg', transparent=True)


fig42, ax42 = new_figure(latex_style=True, font_size=18)
scale = 0.8
fig42.set_size_inches(6.4*scale, 4.8*scale+0.5)

ax42.plot(L_list_sub, M_list_sub2, color='green', linestyle='-', lw=2.5)
ax42.plot(L_list_sup3, M_list_sup3, color='green', linestyle='-', lw=2.5)
ax42.set_xlim(ax22.get_xlim())
ax42.set_ylim(ax22.get_ylim())
ax42.set_yticks([M_list_sup3[-1]])
ax42.grid(zorder=0)
ax42.set_xticks([0,L,2*L])
ax42.set_xticklabels([])
# ax31.set_xticklabels(['Inlet', 'Throat', 'Exit'], color='white')
fig42.tight_layout(pad=0.5)
fig42.savefig('fig42.svg', transparent=True)

fig5, ax5 = new_figure(latex_style=True, box=True)
scale = 0.8
fig5.set_size_inches(6.4*scale, 4.8*scale+0.5)
ax5.plot(L_list_sub, p_p0_list_sub, lw=2.5, color='blue', linestyle='-', label=r'$p/p_{t}$', zorder=4)
ax5.plot(L_list_sup, p_p0_list_sup, lw=2.5, color='blue', linestyle='-', zorder=4)
ax5.set_yticks([round(x,y) for x,y in zip([p_p0_list_sup[-1], p0_p_crit, 1],[3,3,0])])
ax5.set_yticklabels([str(round(x,y)) for x,y in zip([p_p0_list_sup[-1], p0_p_crit, 1],[3,3,0])])
ax5.set_ylabel(r'$p/p_{t}$')
ax5.set_xticks([0,L,2*L])
ax5.set_xticklabels(['Inlet', 'Throat', 'Exit'])
ax5.set_xlim([0,2*L*1.25])
ax5.set_ylim(ax21.get_ylim())
ax5.grid(zorder=0)
fig5.tight_layout(pad=0.5)


fig41, ax41 = new_figure(latex_style=True, font_size=18)
scale = 0.8
fig41.set_size_inches(6.4*scale, 4.8*scale+0.5)
ax41.plot(L_list_sub, p_p0_list_sub2, color='blue', lw=2.5)
ax41.plot(L_list_sup3, p_p0_list_sup3, color='blue', lw=2.5)
# ax41.plot(L_list_sub, p_p0_list_sub, lw=2.5, color='r', linestyle='--', label=r'$p/p_{t}$', zorder=4)
# ax41.plot(L_list_sup, p_p0_list_sup, lw=2.5, color='r', linestyle='--', zorder=4)
ax41.plot([0,2*L],[p_p0_list_sup3[-1], p_p0_list_sup3[-1]], color='black', linestyle='--')
ax41.set_ylim(ax21.get_ylim())
ax41.set_xticks([0,L,2*L])
ax41.set_xticklabels([])
ax41.set_yticks([p_p0_list_sup3[-1]])
ax41.set_yticklabels([np.round(p_p0_list_sup3[-1],3)])
ax41.set_xlim(ax21.get_xlim())
# ax31.set_xticklabels(['Inlet', 'Throat', 'Exit'], color='white')
fig41.tight_layout(pad=0.5)
fig41.savefig('fig41.svg', transparent=True)


print(p_p0_list_sub2)
plt.show()