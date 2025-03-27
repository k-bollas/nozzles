import os
import sys
script_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, script_folder)
from figure_properties import new_figure
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

Mth = 0.75
M2 = 0.3
M_FL = 0.85
P_inf = 23841
T_inf = 218.81
P_t = P_inf*(1+0.2*M_FL**2)**3.5
P_th = P_t*(1+0.2*Mth**2)**(-3.5)
A2_Ath = Mth/M2*((1+0.2*M2**2)/(1+0.2*Mth**2))**(2.4/0.8)

def balance_forces(M0):

    P0 = P_t*(1+0.2*M0**2)**(-3.5)
    A0_Ath = Mth/M0*((1+0.2*M0**2)/(1+0.2*Mth**2))**(2.4/0.8)

    Ainf_A0 = M0/M_FL*((1+0.2*M_FL**2)/(1+0.2*M0**2))**(2.4/0.8)
    Ainf_Ath = Ainf_A0*A0_Ath

    Ath = 1
    A_inf = Ainf_Ath*Ath
    A0 = A0_Ath*Ath

    F_pre = P_inf*A_inf - P_inf*(A_inf-A0) - P0*A0
    F_lip = P0*A0 - P_th*Ath

    print(M0, F_pre, F_lip, F_pre + F_lip)

    return F_pre + F_lip

M0 = brentq(f= balance_forces, a =0.2, b=0.75, xtol=1e-5)
print(M0)