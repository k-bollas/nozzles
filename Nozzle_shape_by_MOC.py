import os
import sys
script_folder = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_folder)
from figure_properties import new_figure
from NozzleMOC__ import NozMOC
import matplotlib.pyplot as plt
import time

if __name__ == '__main__':
    start = time.perf_counter()
    Gamma = 1.4
    Me = 2.5
    n = 25
    y0, x0, y_ref = 1, 0, 1

    [xw, yw, Noz] = NozMOC(Me, n, Gamma, x0, y0, y_ref)

    print(yw[-1])
    x, y, M = Noz['x'], Noz['y'], Noz['M']
    fig, ax = new_figure()

    color1 = 'blue'
    color2 = [68/255, 114/255, 196/255]
    linewidth1 = 4
    linewidth2 = 0.5
    linestyle1 = 'solid'
    linestyle2 = 'solid'

    ax.plot([xw[0]-0.7, xw[0]], [yw[0], yw[0]], color=color1, linewidth=linewidth1, linestyle=linestyle1)
    ax.plot([xw[0]-0.7, xw[0]], [-yw[0], -yw[0]], color=color1, linewidth=linewidth1, linestyle=linestyle1)
    ax.plot(xw, yw, color=color1, linewidth=linewidth1, linestyle=linestyle1)
    ax.plot(xw, -yw, color=color1, linewidth=linewidth1, linestyle=linestyle1)
    ax.plot([xw[0], xw[0]], [-yw[0], yw[0]], color='k', linewidth=1.5, linestyle='--')

    ax.axis([xw[0]-0.7, max(xw), -1.05*yw[-1], 1.05*yw[-1]])
    # ax.set_xlabel("Undimensional X")
    # ax.set_ylabel("Undimensional Y")
    # ax.set_title('Nozzle Shape')

    # Display lines
    if color2 != None:
        for i in range(0,n):
            ax.plot([x0,x[0][i]],[y0,y[0][i]],color=color2, linewidth=linewidth2, linestyle=linestyle2)
            ax.plot([x[i][n-1-i],xw[i+1]],[y[i][n-1-i],yw[i+1]],color=color2, linewidth=linewidth2, linestyle=linestyle2)
            ax.plot([x0,x[0][i]],[-y0,-y[0][i]],color=color2, linewidth=linewidth2, linestyle=linestyle2)
            ax.plot([x[i][n-1-i],xw[i+1]],[-y[i][n-1-i],-yw[i+1]],color=color2, linewidth=linewidth2, linestyle=linestyle2)

        for i in range(0,n-1):
            ax.plot(x[i][0:n-i],y[i][0:n-i],color=color2, linewidth=linewidth2, linestyle=linestyle2)
            ax.plot(x[i][0:n-i],-y[i][0:n-i],color=color2, linewidth=linewidth2, linestyle=linestyle2) 

        for c in range(0,n):
            for r in range(1,n-c):
                ax.plot([x[r][c],x[r-1][c+1]],[y[r][c],y[r-1][c+1]],color=color2, linewidth=linewidth2, linestyle=linestyle2)
                ax.plot([x[r][c],x[r-1][c+1]],[-y[r][c],-y[r-1][c+1]],color=color2, linewidth=linewidth2, linestyle=linestyle2)

    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.subplots_adjust(0.125, 0.11, 0.785, 0.88)
    plt.tight_layout(pad=0.05)
    plt.savefig('nozzle.svg', transparent=True)
    plt.show()