import os
import sys
script_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, script_folder)
from figure_properties import new_figure
import matplotlib.pyplot as plt
import numpy as np


def main():

    gamma = 1.35
    R = 287
    A = 1
    NPR_crit = (1 + (gamma-1)/2)**(gamma/(gamma-1))





    pt6 = 180e3
    Tt6 = 1140

    NPR_r = np.linspace(1, NPR_crit, 200)
    M6 = np.sqrt((NPR_r**((gamma-1)/gamma)-1)*(2/(gamma-1)))
    T6 = Tt6/(1+(gamma-1)/2*M6**2)
    p6  = pt6/NPR_r
    rho6 = p6/(R*T6)
    V6 = M6*np.sqrt(gamma*R*T6)
    m_flow = rho6*A*V6

    fig, ax = new_figure(latex_style=True)
    scale = 0.7
    fig.set_size_inches(6.4*scale, 4.8*scale)

    ax.plot(m_flow, NPR_r, '-b', lw=3)
    ax.plot([np.max(m_flow), np.max(m_flow)], [NPR_crit, 1.5*NPR_crit], '-b', lw=3)
    ax.set_xlabel(r'Mass flow rate, $\dot{m}$')
    ax.set_xticks([0, np.max(m_flow)])
    # ax.axvline(np.max(m_flow), ls='--', lw=1.5, color='k')
    ax.grid()
    ax.set_xticklabels(["0", r'$\dot{m}_{\mathrm{choke}}$'])
    ax.set_ylabel(r'$p_{t6} / p_{\mathrm{amb}}$')
    ax.set_yticks([1, NPR_crit, 1.5*NPR_crit])
    ax.set_yticklabels(['1', r'$\mathrm{NPR}_{\mathrm{crit}}$', r'+$\infty$'])
    mchoke = A*pt6/np.sqrt(Tt6)*np.sqrt(gamma/R)*(1+(gamma-1)/2)**(-(gamma+1)/(2*(gamma-1)))
    print(np.sqrt(gamma/R)*(1+(gamma-1)/2)**(-(gamma+1)/(2*(gamma-1))))
    # ax.scatter(mchoke, 1.2*NPR_crit)
    plt.tight_layout(pad=0.5)

    fig, ax = new_figure(latex_style=True)
    scale = 0.5
    fig.set_size_inches(6.4*scale, 4.8*scale)
    ax.plot(1/NPR_r, m_flow, '-r', lw=3)
    ax.plot([0, 1/NPR_crit], [np.max(m_flow), np.max(m_flow)], '-r', lw=3)
    ax.set_xlim([0,1.1])
    ax.set_ylabel(r'Mass flow rate, $\dot{m}$')
    ax.set_yticks([0, np.max(m_flow)])
    ax.grid()
    ax.set_yticklabels(["0", r'$\dot{m}_{\mathrm{choke}}$'])
    ax.set_xlabel(r'$p_{\mathrm{amb}} / p_{t6}$')
    ax.set_xticks([0, 1/NPR_crit, 1])
    ax.set_xticklabels(['0', r'$\frac{1}{\mathrm{NPR}_{\mathrm{crit}}}$', '1'])
    print(np.round(1/NPR_crit, 3))
    plt.tight_layout(pad=0.5)

    plt.show()

def supersonic_flow():
    pt = 5.73*1e5
    Tt = 1000

    Ae = 1.46
    gamma = 1.4
    R = 287

    pe_pt = np.linspace(1,0.528,100)
    Me = np.sqrt((pe_pt**(-(gamma-1)/gamma)-1)*2/(gamma-1))

    m = Ae*pt/np.sqrt(Tt)*np.sqrt(gamma/R)*Me*(1+(gamma-1)/2*Me**2)**(-(gamma+1)/(2*(gamma-1)))

    plt.rcParams.update({'font.family': 'Arial', 'font.size':14, 'grid.linewidth':1, 'grid.linestyle':'--'})
    fig, ax = plt.subplots()
    ax.plot(Me,m, lw=3, zorder=4)
    ax.grid()
    ax.set_xlabel('Mach number (-)')
    ax.set_ylabel('Mass Flow Rate (kg/s)')
    fig.tight_layout(pad=0.3)
    plt.show()


if __name__ == '__main__':
    main()
    # supersonic_flow()