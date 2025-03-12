from matplotlib.pyplot import rcParams, subplots, show, tight_layout, rc
from numpy import size

def new_figure(*args, **kwargs):
    kwargs_params = {
        'font_size': 16,
        'font_family': 'Times New Roman',
        'latex_style': False,
        'color': 'black',
        'box' : False,
    }
    kwargs_params.update(kwargs)

    # rcParams.update({'font.size': kwargs_params['font_size'], 'font.family': kwargs_params['font_family'], 'legend.loc':'best', 'legend.shadow':True,
    # 'legend.framealpha':1, 'legend.fancybox':False, 'legend.edgecolor':'black', 'grid.linestyle':'--', 'grid.linewidth':1.5,
    # 'grid.alpha':0.5, 'lines.linewidth':kwargs_params['linewidth']})
    if kwargs_params['latex_style'] == True:
        rc('text', usetex=True)
        rc('font', family='serif')
        rcParams.update({'font.size':kwargs_params['font_size'], 'legend.loc':'best', 'legend.shadow':True, 'legend.framealpha':1, 'legend.fancybox':False, 'legend.edgecolor': kwargs_params['color'], 'grid.linestyle':'--', 'grid.linewidth':1.5,
            'grid.alpha':0.5, 'axes.edgecolor':kwargs_params['color'], 'xtick.color':kwargs_params['color'], 'ytick.color':kwargs_params['color'], 'axes.labelcolor': kwargs_params['color'], 'legend.labelcolor': kwargs_params['color'] })
    else:
        rc('font', family=kwargs_params['font_family'])

        rcParams.update({'font.size':kwargs_params['font_size'], 'legend.loc':'best', 'legend.shadow':True, 'legend.framealpha':1, 'legend.fancybox':False, 'legend.edgecolor': kwargs_params['color'],
                         'grid.linestyle':'--', 'grid.linewidth':1.5, 'grid.alpha':0.5, 'axes.edgecolor':kwargs_params['color'], 'xtick.color':kwargs_params['color'], 'ytick.color':kwargs_params['color'], 
                         'axes.labelcolor': kwargs_params['color'], 'legend.labelcolor': kwargs_params['color'] })
        
        rcParams['mathtext.fontset'] = 'custom'
        rcParams['mathtext.rm'] = kwargs_params['font_family']
        rcParams['mathtext.it'] = kwargs_params['font_family']
        rcParams['mathtext.bf'] = kwargs_params['font_family'] 


    fig, ax = subplots(*args)
    fig.patch.set_alpha(0.0)

    if kwargs_params['box'] == False:
        if size(ax) > 1:
            for side in ['right','top']:
                for i in range(0, size(ax)):
                        ax[i].spines[side].set_visible(False)
        else:
            for side in ['right','top']:
                ax.spines[side].set_visible(False)

    return fig, ax

def format_exponent(ax, axis='y'):
    # Change the ticklabel format to scientific format
    ax.ticklabel_format(axis=axis, style='sci', scilimits=(-2, 2))

    # Get the appropriate axis
    if axis == 'y':
        ax_axis = ax.yaxis
        x_pos = 0.0
        y_pos = 1.0
        horizontalalignment='left'
        verticalalignment='bottom'
    else:
        ax_axis = ax.xaxis
        x_pos = 1.0
        y_pos = -0.05
        horizontalalignment='right'
        verticalalignment='top'

    # Run plt.tight_layout() because otherwise the offset text doesn't update
    tight_layout()
    ##### THIS IS A BUG 
    ##### Well, at least it's sub-optimal because you might not
    ##### want to use tight_layout(). If anyone has a better way of 
    ##### ensuring the offset text is updated appropriately
    ##### please comment!

    # Get the offset value
    offset = ax_axis.get_offset_text().get_text()

    if len(offset) > 0:
        # Get that exponent value and change it into latex format
        minus_sign = u'u2212'
        expo = float(offset.split('e')[-1][-1])
        offset_text = r'$\times$ 10$^{-%d}$' %expo

        # Turn off the offset text that's calculated automatically
        ax_axis.offsetText.set_visible(False)

        # Add in a text box at the top of the y axis
        ax.text(x_pos, y_pos, offset_text, transform=ax.transAxes,
               horizontalalignment=horizontalalignment,
               verticalalignment=verticalalignment)
    return ax

if __name__ == '__main__':
    fig, ax = new_figure()