import WrightTools as wt
import numpy as np
import matplotlib.pyplot as plt

from . import fit as fit


def PL_g2_plot(ax, d):
    ax.plot(d, linewidth=1)
    ax.fill_between(d.delay.points, d.counts.points, alpha=.3)
    ax.set_xlim(d.delay.min(), d.delay.max())
    ax.set_ylim(0, d.counts.max()*1.05)
    text = 'Area ratio: ' + "{:.2e}".format(d.attrs['arearatio'])
    bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
    ax.text(0.99, .96, text, horizontalalignment='right',
            verticalalignment='top', transform=ax.transAxes, bbox=bbox)


def PL_picotime_plot(ax, d, fitting=True):
    x = d.delay.points
    y = d.counts.points
    ax.plot(d, linewidth=1, alpha=1)
    ax.fill_between(x, y, alpha=.3)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min()+1, y.max()*1.1)
    ax.set_yscale('log')
    if fitting:
        pfit, perr, ymodel = fit.exp_fit(x, y)
        ax.plot(x,ymodel, color='C1', linewidth=4, alpha=.75)
        t1 = pfit[0]  
        text = '$\\mathsf{\\tau_1 =' + str(round(np.abs(t1),1)) + '\\;ns}$'
        bbox = dict(boxstyle='round', fc='C1', ec='C1', alpha=0.5)
        ax.text(0.99, .96, text, horizontalalignment='right',
                verticalalignment='top', transform=ax.transAxes, bbox=bbox)     
        pfit, perr, ymodel = fit.biexp_fit(x, y)
        ax.plot(x,ymodel, color='C2', linewidth=4, alpha=.75)
        text = '$\\mathsf{\\tau_1 =' + str(round(np.abs(pfit[0]),1)) + '\\;ns, \\; \\tau_2 =' + str(round(np.abs(pfit[2]),1)) +'\\;ns, \\;'  + 'A_2/A_1 =' + str(round(np.abs(pfit[3]/pfit[1]),2))+'}$'
        bbox = dict(boxstyle='round', fc='C2', ec='C12', alpha=0.5)
        ax.text(0.99, .8, text, horizontalalignment='right',
                verticalalignment='top', transform=ax.transAxes, bbox=bbox) 
      

def PL_macrotime_plot(ax, d, col):
    ax.plot(d, linewidth=1)
    ax.fill_between(d.labtime.points, d.counts.points, alpha=.3)
    ax.set_xlim(d.labtime.min(), d.labtime.max())
    ax.set_ylim(d.counts.min()+1, d.counts.max()*1.05)
    text = 'Records: ' + "{:.2e}".format(col.attrs['Records'])
    bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
    ax.text(0.99, .96, text, horizontalalignment='right',
            verticalalignment='top', transform=ax.transAxes, bbox=bbox)
    

def PL_fig_plot(col, fitting=True):
    fig, gs = wt.artists.create_figure(width='double', nrows=3, default_aspect=.25, hspace=.7)
    axs = [plt.subplot(gs[i]) for i in range(3)]
    ylabels = ['$\\mathsf{counts \\; [Hz]}$', '$\\mathsf{counts}$', '$\\mathsf{cross-correlation}$']
    xlabels = ['$\\mathsf{labtime \\; [s]}$', '$\\mathsf{delay \\; [ns]}$', '$\\mathsf{\\tau/\\tau_{rep}}$' ]
    for ax, xlabel, ylabel in zip(axs, xlabels, ylabels):
        ax.grid()
        wt.artists.set_ax_labels(ax=ax, xlabel=xlabel, ylabel=ylabel)
    PL_macrotime_plot(axs[0], col.macrohist, col)
    PL_picotime_plot(axs[1], col.picohist, fitting)
    PL_g2_plot(axs[2], col.g2hist)
    axs[0].set_title(col.attrs['identifier'])
    return fig, gs
