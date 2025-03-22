import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def cpr_line():
    AR1 = 1.75
    AR2 = 3.5
    CPR1 = 0.6
    CPR2 = 0.8
    a_ = (CPR2-CPR1)/(AR2-AR1)
    b_ = CPR2 - a_*AR2
    return a_, b_

def lambda_function():
    AR_ = np.linspace(1.5, 3.5, 200)
    a, b_ =cpr_line()
    CPR = a_*AR_ + b_
    CPRI = 1 - (1/AR)**2
    l_factor = CPRI - CPR
    l_func = np.polyfit(AR_, l_factor, 4)
    return l_func


M1 = 0.3
M2 = 0.15
T1 = 288
p1 = 101325
V1 = M1*np.sqrt(1.4*287*T1)
rho1 = p1/(287*T1)
q1 = 1/2*rho1*V1**2
p01	= p1 + q1

p_loss = 0
error = 1
while error > 1e-3:

    p02 = p01 - p_loss
    p2 = p02/(1+1.4/2*M2**2)
    T2 = p2/(rho1*287)
    V2 = M2*np.sqrt(1.4*287*T2)
    C_PR = (p2-p1)/q1
    a_, b_ = cpr_line()
    AR = (C_PR-b_)/a_
    l_factor = np.polyval(lambda_function(), AR)
    p_loss_new = l_factor*q1
    error = np.abs(p_loss_new - p_loss)
    p_loss = p_loss_new

    print(l_factor, error)