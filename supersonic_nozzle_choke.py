import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sc
from figure_properties import new_figure
from scipy.optimize import fsolve

# ----------------------
A_star = 1
p_total = 5.73e5
T_total = 1400
A_e_Astar = 1.44
A_e = A_e_Astar*A_star

p_e_st = 101325

# ----------------------

pe_p_ratio = np.linspace(1, p_e_st/p_total, 100)
Me, m_flow, M_th = [np.zeros(len(pe_p_ratio)) for _ in range(3)]

for i, _ in enumerate(pe_p_ratio):
    Me[i] = np.sqrt(2/0.4*(pe_p_ratio[i]**(-1/3.5)-1))
    m_flow[i] = A_e*p_total/np.sqrt(T_total)*np.sqrt(1.4/287)*Me[i]*(1+0.2*Me[i]**2)**(-2.4/0.8)
    M_th[i], = fsolve(func= lambda x: A_star*p_total/np.sqrt(T_total)*np.sqrt(1.4/287)*x*(1+0.2*x**2)**(-2.4/0.8) - m_flow[i], x0=0.1, xtol=1e-5)
    if np.round(M_th[i],2) - 1 == 0: M_th[i] = np.round(M_th[i],2)


    # if M_th[i] > 1 : 
    #     M_th[i] = 1
    #     m_flow[i] = A_star*p_total/np.sqrt(T_total)*np.sqrt(1.4/287)*1*(1+0.2)**(-2.4/0.8)


print(Me, m_flow, M_th)