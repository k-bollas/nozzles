import os
import sys
script_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, script_folder)
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import colorsys
import numpy as np
import time

import numpy as np

def PMF(M,Gamma):
    Gp=Gamma+1 # γ-1
    Gm=Gamma-1 # γ+1
    if Gamma>1:
        nu = np.sqrt(Gp/Gm)*np.arctan(np.sqrt(Gm*(np.power(M,2)-1)/Gp))-np.arctan(np.sqrt(np.power(M,2)-1)) # Pranldt-Meyer angle in rad
        nu = np.rad2deg(nu) # Pranldt-Meyer angle in degrees
    elif Gamma==1:
        nu = np.arctan(1/(np.sqrt(np.power(M,2)-1)))+np.deg2rad(np.sqrt(np.power(M,2)-1))-np.deg2rad(45) # Pranldt-Meyer angle in rad
        nu = np.rad2deg(nu) # Pranldt-Meyer angle in degrees
    else:
        nu = complex(0,-1)/np.sqrt(complex(Gm/Gp,0))*np.arctanh(-np.sqrt(complex(Gm/Gp,0))/complex(0,1)*np.sqrt(np.power(M,2)-1))-np.arctan(np.sqrt(np.power(M,2)-1)) # Pranldt-Meyer angle in rad
        assert(np.imag(nu.any())==0),"Prandlt-Mayer nu angle calculation error : The imaginary part is not zero" # Stop excecution if Prandlt-Meyer angle is not real number
        nu = np.rad2deg(np.real(nu)) # Pranldt-Meyer angle in degrees
    mu = np.arcsin(1/M) # Mach angle in rad
    mu = np.rad2deg(mu) # Mach angle in degrees
    return M,nu,mu

def goal_seek(target,_threshold,Gamma): # Gold seek algorithm for finding the proper Mach number for a know Prandlt-Meyer angle
    Mmax = 4; Mmin = 1
    numax = PMF(Mmax,Gamma)[1]; numin = PMF(Mmin,Gamma)[1]
    fmax = target - numax; fmin = target - numin
    Mprev = Mmin; fprev = fmin; 
    L_ = fmin; R_ = fmax; Range = abs(Mmax-Mmin)
    while Range >= _threshold:
        Mnext = (Mmin*R_-Mmax*L_)/(R_-L_)
        nunext = PMF(Mnext,Gamma)[1]; fnext = target - nunext
        if fnext == 0:
            break
        elif fmin*fnext<=0:
            Mmax = Mnext; fmax = fnext; R_ = fnext
            if fprev*fnext>0:
                L_ = L_/2
        else:
            Mmin = Mnext; fmin = fnext; L_ = fnext
            if fprev*fnext>0:
                R_ = R_/2
        Mprev = Mnext; fprev = fnext; Range = (abs(Mmax-Mmin))
#         print(f'Range is: {Mmin}  ----  {Mnext}   ----   {Mmax}, ----> {fnext} --------------------------------',end="\r", flush=True)
    return Mnext

def secant(x1,x2,f,tol, *args):
    err = 1
    stepNum = 0
    while err > tol:
        stepNum=stepNum+1
        A = np.array([[1,1],[f(x1), f(x2)]])
        B = np.array([[1], [0]])
        p = np.matmul(np.linalg.inv(A), B)
        x_n = float(sum([k1*k2 for k1,k2 in zip(p, [x1,x2])]))
        if 'show' in args:
            print(x_n, x1, x2)
        if x_n < 0 and 'only_possitive' in args:
            x_n = float(np.random.rand(1)*(x1 + x2)/2)
        err = abs(x2 - x_n)
        x1,x2 = x2, x_n
        # print(x1, x2)
    
    return x_n, f(x_n), stepNum


# Inputs
def MNM(M,nu,mu,Gamma): # Mach-Nu-Mu Function: One input only must be a number, the other values must be zero
    # ex: M = 2, mu = 0, nu = 0, Gamma = 1.4 or M = 0,mu = 30, nu = 0, Gamma 1.4 etc.
    M = np.array(M)
    mu = np.array(mu)
    nu = np.array(nu)

    Error = True
    if M.any()!=0 and nu.any()==0 and mu.any()==0:
        Error = False
        [M,nu,mu] = PMF(M,Gamma)
    elif mu.any()!=0 and M.any()==0 and nu.any()==0:
        Error = False
        M = 1/np.sin(np.deg2rad(mu))
        [M,nu,mu] = PMF(M,Gamma)
    elif nu.any()!=0 and M.any()==0 and mu.any()==0:
        Error = False
        nu = [nu] if np.size(nu) == 1 else nu
        M = np.array([secant(1, 1.5, lambda x: PMF(x,Gamma)[1] - y, 1e-5, 'only_possitive')[0] for y in nu])
        mu = np.rad2deg(np.arcsin(1/M))    
        if np.size(M) == 1:
            M, nu, mu = M[0], nu[0], mu[0]

        # if np.size(nu)==1:
        #     M = goal_seek(nu,1e-6,Gamma)
        #     mu = np.arcsin(1/M)
        #     mu = np.rad2deg(mu)
        # else:
        #     M = np.ones(len(nu))
        #     mu = np.ones(len(nu))
        #     for i in range(0,len(nu)):
        #         M[i] = goal_seek(nu[i],1e-6,Gamma)
        #         mu[i] = np.arcsin(1/M[i])
        #         mu[i] = np.rad2deg(mu[i])
    else:
        if Error == True:
            raise ValueError("Error! Only one input must be not zero - the other 2 must be equal with zero")
    return M,nu,mu


def NozMOC(Me, n, Gamma, x0, y0, y_ref, *args) -> list:
    # Me ---> Exit Mach number
    # n --> Number of Characteristic lines
    # Gamma --> Ratio of heat capacities γ = Cp/Cv
    # x0 , y0 --> Coordinates at the throat
    # y_ref --> Reference coordinates for dimensionless coordinates

    if Me < 1:
        raise ValueError(' Subsonic Flow. Exit Mach number must be above 1 (>1)')
    if int(n) < 1:
        raise ValueError(' We have to set at least one characteristic. Number of characteristics should be not less than 1 (<1)')
    if n != int(n):
        raise ValueError(' Number of characteristics must be an integer number')

    B = MNM(Me,0,0,Gamma)[1]
    ThetaMax = B/2
    dT = ThetaMax/n
    
    # Throat coordinates
    y0, x0 = y0/y_ref, x0/y_ref

    Theta = np.zeros([n,n])
    Theta[0] = np.linspace(dT,ThetaMax,n)
    Nu = np.zeros([n,n])
    Nu[0] = Theta[0]
    Km = np.zeros([n,n])
    Km[0] = Nu[0]+Theta[0]
    Kp = np.zeros([n,n])
    Kp[0] = Nu[0]-Theta[0]
    M = np.zeros([n,n])
    M[0] = MNM(0,Nu[0],0,Gamma)[0]
    Mu = np.zeros([n,n])
    Mu[0] = MNM(0,Nu[0],0,Gamma)[2]

    # Characterictic lines
    x, y = np.zeros([n,n]), np.zeros([n,n])
    # First characterictic line
    y[0][0] = 0
    x[0][0] = x0-y0/np.tan(np.deg2rad(Theta[0][0]-Mu[0][0]))
    for i in range(1,n):
        s1 = np.tan(np.deg2rad(Theta[0][i]-Mu[0][i])) # Slope of C- line
        s2 = np.tan(np.deg2rad((Theta[0][i-1]+Mu[0][i-1]+Theta[0][i]+Mu[0][i])/2)) # Slope of C+ line
        x[0][i]=((y[0][i-1]-x[0][i-1]*s2)-(y0-x0*s1))/(s1-s2) # Intersection point of the two characterictis
        y[0][i]=y[0][i-1]+(x[0][i]-x[0][i-1])*s2 # Intersection point of the two characterictis

    # Rest of the characterictic lines
    for j in range(1,n):
        if 'display_progress' in args: print("{} % completed.....".format(round(j/n*100,2)), end="\r", flush=True)
        for i in range(0,1+n-(j+1)):
            Km[j][i] = Km[j-1][i+1]
            if i == 0:
                Theta[j][i] = 0
                Kp[j][i] = -Km[j][i]
                Nu[j][i] = Km[j][i]
                [M[j][i],Nu[j][i],Mu[j][i]] = MNM(0,Nu[j][i],0,Gamma)
                s1 = np.tan(np.deg2rad((Theta[j-1][i+1]-Mu[j-1][i+1]+Theta[j][i]-Mu[j][i])/2))
                x[j][i] = x[j-1][i+1] - y[j-1][i+1]/s1
                y[j][i] = 0
            else:
                Kp[j][i] = Kp[j][i-1]
                Theta[j][i] = (Km[j][i]+Kp[j][i])/2
                Nu[j][i] = (Km[j][i]-Kp[j][i])/2
                [M[j][i],Nu[j][i],Mu[j][i]] = MNM(0,Nu[j][i],0,Gamma)
                s1 = np.tan(np.deg2rad((Theta[j-1][i+1]-Mu[j-1][i+1]+Theta[j][i]-Mu[j][i])/2))
                s2 = np.tan(np.deg2rad((Theta[j][i-1]+Mu[j][i-1]+Theta[j][i]+Mu[j][i])/2))
                x[j][i] = ((y[j][i-1]-x[j][i-1]*s2)-(y[j-1][i+1]-x[j-1][i+1]*s1))/(s1-s2)
                y[j][i] = y[j][i-1]+(x[j][i]-x[j][i-1])*s2

    # Coordinates of the wall of the nozzle
    if 'display_progress' in args: print("{} % completed.....".format(round(100,2)), end="\r", flush=True)
    xw = np.zeros(n+1); yw = np.zeros(n+1)

    #Throat coordinates
    xw[0], yw[0] = x0, y0
    s11 = np.tan(np.deg2rad(ThetaMax))
    s22 = np.tan(np.deg2rad(Theta[0][n-1]+Mu[0][n-1]))
    # First characterictic coordinate points
    xw[1] = ((y[0][n-1]-x[0][n-1]*s22)-(yw[0]-xw[0]*s11))/(s11-s22)
    yw[1] = yw[0]+(xw[1]-xw[0])*s11
    # Rest coordinate points
    for j in range(2,n+1):
        s11 = np.tan(np.deg2rad((Theta[j-2][n-j+1]+Theta[j-1][n-j])/2))
        s22 = np.tan(np.deg2rad(Theta[j-1][n-j]+Mu[j-1][n-j]))
        xw[j] = ((y[j-1][n-j]-x[j-1][n-j]*s22)-(yw[j-1]-xw[j-1]*s11))/(s11-s22)
        yw[j] = yw[j-1]+((xw[j]-xw[j-1])*s11)

    Noz = {
        'x': x,
        'y': y,
        'M': M,
        'Theta': Theta,
        'Nu': Nu,
        'Km': Km,
        'Kp': Kp,
        'Mu': Mu,
    }
    return xw,yw,Noz

def Plot_of_Nozzle(xw, yw, Noz, **kwargs):
    # Set Plot Properties :
    # ----- color1: --> matplotlib.color either 'string' or [R/255, G/255, B/255] --> Color of the nozzle wall
    # ----- color2: matplotlib.color either 'string' or [R/255, G/255, B/255] --> Color of the nozzle characteristic lines
    # ----- linewidth1: --> float --> Line Width of the nozzle wall
    # ----- linewidth2: --> float --> Line Width of the nozzle characteristic lines
    # ----- linestyle1: --> 'string' --> Line Style of the nozzle wall
    # ----- linestyle2: --> 'string' --> Line Style of the nozzle characteristic lines
    # ----- symmetry_axis: --> True or False --> Show or hide symmetry axis in the plot
    # ----- display: --> 'string' --> Insert value 'half' if you want to display half of the nozzle plot

    color1 = kwargs.get('color1') if 'color1' in kwargs else 'blue'
    color2 = kwargs.get('color2') if 'color2' in kwargs else [0.26, 0.447, 0.768]
    font_size = kwargs.get('font_size') if 'font_size' in kwargs else 16
    font_family = kwargs.get('font_family') if 'font_family' in kwargs else 'Times New Roman'
    linewidth1 = kwargs.get('linewidth1') if 'linewidth1' in kwargs else 3
    linewidth2 = kwargs.get('linewidth2') if 'linewidth2' in kwargs else 0.5
    linestyle1 = kwargs.get('linestyle1') if 'linestyle1' in kwargs else 'solid'
    linestyle2 = kwargs.get('linestyle2') if 'linestyle2' in kwargs else 'solid'
    symmetry_axis = kwargs.get('symmetry_axis') if 'symmetry_axis' in kwargs else True
    display = kwargs.get('display') if 'display' in kwargs else 'whole'

    x, y, M = Noz['x'], Noz['y'], Noz['M']
    plt.rcParams['font.family'] = font_family
    plt.rcParams['font.size'] = font_size

    fig, ax = plt.subplots()
    ax.plot(xw, yw, color=color1, linewidth=linewidth1, linestyle=linestyle1)
    if not(display) == 'half': 
        ax.plot(xw, -yw, color=color1, linewidth=linewidth1, linestyle=linestyle1)
        ax.axis([xw[0], max(xw), -1.05*yw[-1], 1.05*yw[-1]])
    else:
        ax.axis([xw[0], max(xw), 0, 1.05*yw[-1]])
    if symmetry_axis==True: ax.plot([xw[0], max(xw)], [0,0], color='black', linewidth=1, linestyle='--')
    ax.set_xlabel("Undimensional X")
    ax.set_ylabel("Undimensional Y")
    ax.set_title('Nozzle Shape')

    # Display lines
    if color2 != None:
        for i in range(0,n):
            ax.plot([x0,x[0][i]],[y0,y[0][i]],color=color2, linewidth=linewidth2, linestyle=linestyle2)
            ax.plot([x[i][n-1-i],xw[i+1]],[y[i][n-1-i],yw[i+1]],color=color2, linewidth=linewidth2, linestyle=linestyle2)
            if not(display) == 'half':
                ax.plot([x0,x[0][i]],[-y0,-y[0][i]],color=color2, linewidth=linewidth2, linestyle=linestyle2)
                ax.plot([x[i][n-1-i],xw[i+1]],[-y[i][n-1-i],-yw[i+1]],color=color2, linewidth=linewidth2, linestyle=linestyle2)

        for i in range(0,n-1):
            ax.plot(x[i][0:n-i],y[i][0:n-i],color=color2, linewidth=linewidth2, linestyle=linestyle2)
            if not(display) == 'half': ax.plot(x[i][0:n-i],-y[i][0:n-i],color=color2, linewidth=linewidth2, linestyle=linestyle2) 

        for c in range(0,n):
            for r in range(1,n-c):
                ax.plot([x[r][c],x[r-1][c+1]],[y[r][c],y[r-1][c+1]],color=color2, linewidth=linewidth2, linestyle=linestyle2)
                if not(display) == 'half': ax.plot([x[r][c],x[r-1][c+1]],[-y[r][c],-y[r-1][c+1]],color=color2, linewidth=linewidth2, linestyle=linestyle2)

    # ax.set_yticks([])
    # ax.set_xticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    plt.subplots_adjust(0.125, 0.11, 0.785, 0.88)

    return fig, ax

def Mach_Contour_of_Nozzle(xw, yw, Noz, **kwargs):
    # Set Plot Properties :
    # ----- levels: --> float --> Select the levels of contours in the plot
    # ----- colormap: --> matplotlib.colormap 'string' --> Type of the colormap
    # ----- colorbar: --> True or False --> Show or Hide colorbar in the plot
    # ----- scatter: --> True or False --> Show or Hide calculation points.
    # ----- color1: --> matplotlib.color either 'string' or [R/255, G/255, B/255] --> Color of the nozzle wall
    # ----- color2: matplotlib.color either 'string' or [R/255, G/255, B/255] --> Color of the nozzle characteristic lines
    # ----- linewidth1: --> float --> Line Width of the nozzle wall
    # ----- linewidth2: --> float --> Line Width of the nozzle characteristic lines
    # ----- linestyle1: --> 'string' --> Line Style of the nozzle wall
    # ----- linestyle2: --> 'string' --> Line Style of the nozzle characteristic lines
    # ----- symmetry_axis: --> True or False --> Show or hide symmetry axis in the plot
    # ----- display: --> 'string' --> Insert value 'half' if you want to display half of the nozzle plot

    def man_cmap(cmap, value=1.):
        colors = cmap(np.arange(cmap.N))
        hls = np.array([colorsys.rgb_to_hls(*c) for c in colors[:,:3]])
        hls[:,1] *= value
        rgb = np.clip(np.array([colorsys.hls_to_rgb(*c) for c in hls]), 0,1)
        return mcolors.LinearSegmentedColormap.from_list("", rgb)


    color1 = kwargs.get('color1') if 'color1' in kwargs else 'black'
    color2 = kwargs.get('color2') if 'color2' in kwargs else 'black'
    font_size = kwargs.get('font_size') if 'font_size' in kwargs else 14
    font_family = kwargs.get('font_family') if 'font_family' in kwargs else 'Times New Roman'
    linewidth1 = kwargs.get('linewidth1') if 'linewidth1' in kwargs else 3
    linewidth2 = kwargs.get('linewidth2') if 'linewidth2' in kwargs else 0.2
    linestyle1 = kwargs.get('linestyle1') if 'linestyle1' in kwargs else 'solid'
    linestyle2 = kwargs.get('linestyle2') if 'linestyle2' in kwargs else 'dashed'
    symmetry_axis = kwargs.get('symmetry_axis') if 'symmetry_axis' in kwargs else True
    display = kwargs.get('display') if 'display' in kwargs else 'whole'
    colormap = kwargs.get('colormap') if 'colormap' in kwargs else 'jet'
    colorbar = kwargs.get('colorbar') if 'colorbar' in kwargs else True
    levels = kwargs.get('levels') if 'levels' in kwargs else 20
    scatter = kwargs.get('scatter') if 'scatter' in kwargs else False

    fig, ax = Plot_of_Nozzle(xw, yw, Noz, color1=color1, color2=color2, font_size=font_size, font_family=font_family, 
    linewidth1=linewidth1, linewidth2=linewidth2, linestyle1=linestyle1, linestyle2=linestyle2, symmetry_axis=symmetry_axis, display=display)
    ax.set_title('Nozzle Mach Number Contour')

    x, y, M = Noz['x'], Noz['y'], Noz['M']
    Theta, Mu = Noz['Theta'], Noz['Mu']
    x_c, y_c, M_c = [], [], []    
    x_c, y_c, M_c = np.zeros([len(x[0])+2,len(x[0])+2]), np.zeros([len(x[0])+2,len(x[0])+2]), np.zeros([len(x[0])+2,len(x[0])+2])
    x_c[0] = np.linspace(0, x0, len(x[0])+2)
    y_c[0] = np.linspace(0, y0, len(x[0])+2)
    M_c[0] = np.ones(len(x[0])+2)

    x_c[1] = np.ones(len(x[0])+2)*x0 + 1e-6*(max(xw)-x0)
    y_c[1] = np.zeros(len(x[0])+2)
    y_c[1][1:-1] = np.array([np.tan(np.deg2rad(Theta[0][i-1] - Mu[0][i-1]))*(x_c[1][i] - x0) + y0 for i in range(1,len(x[0])+1)])
    y_c[1][-1] = np.polyfit([xw[0], xw[1]],[yw[0], yw[1]], 1)[0]*(x_c[1][-1] - x0) + y0
    M_c[1][0] = 1
    M_c[1][1:-1] = M[0]
    M_c[1][-1] = M[0][-1]

    for i in range(0, len(x[0])-1):
        x_c[i+2][0] = x[i][0]
        y_c[i+2][0] = y[i][0]
        M_c[i+2][0] = M[i][0]
        x_c[i+2][1:-(i+2)] = x[i][0:-(i+1)]
        y_c[i+2][1:-(i+2)] = y[i][0:-(i+1)]
        M_c[i+2][1:-(i+2)] = M[i][0:-(i+1)]
        x_c[i+2][-(i+2):] = np.linspace(x[i][-(i+1)], xw[(i+1)], i+2)
        y_c[i+2][-(i+2):] = np.linspace(y[i][-(i+1)], yw[(i+1)], i+2)
        M_c[i+2][-(i+2):] = np.ones(i+2)*M[i][-(i+1)]

    x_c[-1] = np.linspace(x[-1][0], xw[-1], len(x[0])+2)
    y_c[-1] = np.linspace(y[-1][0], yw[-1], len(x[0])+2)
    M_c[-1] = np.ones(len(x[0])+2)*M[-1][0]

    x_c, y_c, M_c = list(x_c), list(y_c), list(M_c)
    x_c.append(np.linspace(xw[-1], xw[-1], len(x[0])+2))
    y_c.append(np.linspace(0, yw[-1], len(x[0])+2))
    M_c.append(np.ones(len(x[0])+2)*M[-1][0]+1e-6)

    x_c, y_c, M_c = np.array(x_c), np.array(y_c), np.array(M_c)

    cax = ax.contourf(x_c, y_c, M_c, levels, cmap=man_cmap(plt.get_cmap(colormap), 1.0))
    if not(display) == 'half': ax.contourf(x_c, -y_c, M_c, levels, cmap=man_cmap(plt.get_cmap(colormap), 1.0))
    if scatter == True:
        ax.scatter(x_c,y_c, color='black')
        if not(display) == 'half': ax.scatter(x_c,-y_c, color='black')
    if colorbar == True: 
        c1 = fig.colorbar(cax, ax=ax, ticks=np.linspace(1,round(np.amax(M),2),11))
        c1.set_label('Mach Number')


    plt.subplots_adjust(0.125, 0.11, 0.95, 0.88)

    return fig, ax

if __name__ == '__main__':
    start = time.perf_counter()
    Gamma = 1.4
    Me = 2.5
    n = 20
    y0, x0, y_ref = 1, 0, 1

    [xw, yw, Noz] = NozMOC(Me, n, Gamma, x0, y0, y_ref)

    print(yw[-1])
    x, y, M = Noz['x'], Noz['y'], Noz['M']
    fig1, ax1 = Plot_of_Nozzle(xw, yw, Noz)
    fig2, ax2 = Mach_Contour_of_Nozzle(xw, yw, Noz, levels=100, color1='#e65c00', colormap='inferno')
    end = time.perf_counter()
    print(' ---- Time Elapsed : {:.3f} seconds '.format(end-start))
    plt.show()