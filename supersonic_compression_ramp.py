import os
import sys
script_folder = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_folder)
from figure_properties import new_figure
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from ambiance.ambiance import Atmosphere

def mass_flow_rate_example():
    M = 2.0
    atm = Atmosphere(18290)
    R=287
    gamma=1.4
    P = float(atm.pressure)
    T = float(atm.temperature)
    print(P/1e3, T)
    rho = P/(R*T)
    c = np.sqrt(gamma*R*T)
    V = M*c
    A = 2.5
    W = np.sqrt(A)

    m_flow = rho*A*V
    print(W, m_flow)



mass_flow_rate_example()

sys.exit()

M1 = np.linspace(1.0, 4.0)
gamma = 1.4
Cp = 1005
R = 287


def normal_shock(M1):

    M2 = np.sqrt((1+(gamma-1)/2*M1**2)/(gamma*M1**2-(gamma-1)/2))
    p2p1 = 1 + 2*gamma/(gamma+1)*(M1**2-1)
    T2T1 = p2p1*(2+(gamma-1)*M1**2)/((gamma+1)*M1**2)
    Ds = Cp*np.log(T2T1) - R*np.log(p2p1)
    pt2pt1 = np.exp(-Ds/R)
    return M2, p2p1, T2T1, Ds, pt2pt1

def theta_beta_M(b,theta,M):
    if b > 90: b = 1.0
    return np.tan(np.deg2rad(theta)) - 2/np.tan(np.deg2rad(b))*(M**2*(np.sin(np.deg2rad(b)))**2-1)/(M**2*(gamma+np.cos(np.deg2rad(2*b)))+2)

def compressor_ramps(n):

    M1_ = np.linspace(1.0, 4.0) if n == 1 else np.linspace(1.0, 4.0)

    M, Mn, p2p1, T2T1, Ds, pt2pt1 = [np.zeros([n+1, np.size(M1_)]) for _ in range(6)]
    beta = np.zeros([n, np.size(M1_)])
    M[0,:] = M1_
    if n > 1:
        theta1 = 12 # deg
        theta = np.arange(theta1/(n-1), theta1+0.01, theta1/(n-1))
        for i in range(0,n-1):
            beta[i,:] = np.array([fsolve(lambda x: theta_beta_M(x,theta[i],m), x0=45, xtol=1e-5)[0] for m in M[i,:]])
            Mn[i,:] = M[i,:]*np.sin(np.deg2rad(beta[i,:]))
            Mn[i+1,:], p2p1[i+1,:], T2T1[i+1,:], Ds[i+1,:], pt2pt1[i+1,:] = normal_shock(Mn[i,:])
            M[i+1,:] = Mn[i+1,:]/np.sin(np.deg2rad(beta[i,:]-theta[i]))
            print(beta[i,0], M[i,0], Mn[i,0], Mn[i+1,0], M[i+1,0])
    M[n,:], p2p1[n,:], T2T1[n,:], Ds[n,:], pt2pt1[n,:] = normal_shock(M[n-1,:])

    pt2pt1 = np.delete(pt2pt1, 0, axis=0).squeeze()
    p2p1 = np.delete(p2p1, 0, axis=0).squeeze()
    T2T1 = np.delete(T2T1, 0, axis=0).squeeze()
    Ds = np.delete(Ds, 0, axis=0).squeeze()

    if n > 1:
        pt1pt0 = np.prod(pt2pt1, axis=0)
    else:
        pt1pt0 = pt2pt1

    return M, p2p1, T2T1, Ds, pt1pt0

# Case 1 - Just normal Shock


Case1 = {}
Case1['M'], _, _, _, Case1['pt1pt0'] = compressor_ramps(1)

# Case 2 - Normal Shock + 1 Oblque shock
Case2 = {}
Case2['M'], _, _, _, Case2['pt1pt0'] = compressor_ramps(2)

# Case 3 - Normal Shock + 2 Oblque shocks
Case3 = {}
Case3['M'], _, _, _, Case3['pt1pt0'] = compressor_ramps(3)

# Case 4 - Normal Shock + 4 Oblque shocks
Case4 = {}
Case4['M'], _, _, _, Case4['pt1pt0'] = compressor_ramps(4)

# Case 5
Case5 = {}
Case5['M'], _, _, _, Case5['pt1pt0'] = compressor_ramps(10)


# Case5 = {}
# Case5['M'], _, _, _, Case5['pt1pt0'] = compressor_ramps(200)

fig, ax = new_figure()
scale = 0.8
fig.set_size_inches(6.4*scale, 4.8*scale)

ax.plot(Case1['M'][0,:],  Case1['pt1pt0'], '-b', lw=2.5, zorder=3)
ax.plot(Case2['M'][0,:],  Case2['pt1pt0'], '-b', lw=2.5, zorder=3)
ax.plot(Case3['M'][0,:],  Case3['pt1pt0'], '-b', lw=2.5, zorder=3)
ax.plot(Case4['M'][0,:],  Case4['pt1pt0'], '-b', lw=2.5, zorder=3)
ax.plot(Case5['M'][0,:],  Case5['pt1pt0'], '-b', lw=2.5, zorder=3)

ax.set_xlabel('Flight Mach number')
ax.set_ylabel(r'$p_{t1}/p_{t0}$')
ax.grid(lw=1, ls=':', zorder=0)
fig.tight_layout(pad=0.5)
plt.show()