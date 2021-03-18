import WrightTools as wt
import numpy as np
import matplotlib.pyplot as plt

from . import fit as fit


def PL_g2_plot(ax, d):
    ax.plot(d, linewidth=1)
    ax.fill_between(d.delay.points, d.counts.points, alpha=0.3)
    ax.set_xlim(d.delay.min(), d.delay.max())
    # ax.set_ylim(0, d.counts.max() * 1.05)
    text = "Area ratio: " + "{:.2e}".format(d.attrs["arearatio"])
    bbox = dict(boxstyle="round", fc="blanchedalmond", ec="orange", alpha=0.5)
    ax.text(
        0.99,
        0.96,
        text,
        horizontalalignment="right",
        verticalalignment="top",
        transform=ax.transAxes,
        bbox=bbox,
    )


def PL_picotime_plot(ax, d, fitting=1):
    maxes = 10
    mins = 1e4
    for chan, color, i in zip(d.channels, ["C0", "C1", "C2"], [0, 1, 2]):
        x = d.delay.points
        y = chan.points
        if y.max() > maxes:
            maxes = y.max()
        if y.min() < mins:
            mins = y.min()
        ax.plot(x, y, linewidth=1, alpha=.5, color=color)
        # ax.fill_between(x, y, alpha=0.3, color=color)
        if fitting in [1, 2, 3]:
            ts = []
            if fitting == 1:
                pfit, perr, ymodel = fit.exp_fit(x, y)
                t1, A1, B, t0, fwhm = pfit
                ts.append(t1)
            if fitting == 2:
                pfit, perr, ymodel = fit.biexp_fit(x, y)
                t1, A1, t2, A2, B, t0, fwhm = pfit
                ts.append(t1)
                ts.append(t2)
            if fitting == 3:
                pfit, perr, ymodel = fit.triexp_fit(x, y)
                t1, A1, t2, A2, t3, A3, B, t0, fwhm = pfit
                ts.append(t1)
                ts.append(t2)
                ts.append(t3)
            label = "$\\mathsf{\\tau_i ="
            for t in ts:
                label = label + "\\;" + str(round(np.abs(t), 1)) + ","
            label = label[:-1] + "\\;ns}$"
            ax.plot(x, ymodel, color=color, linewidth=3, alpha=1, label=label)
            ax.legend()
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(mins + 1, maxes * 1.1)
    ax.set_yscale("log")
    


def PL_macrotime_plot(ax, d, col):
    for chan in d.channels:
        ax.plot(d.labtime.points, chan.points, linewidth=2)
        # ax.set_ylim(chan.points.min() + 1, chan.points.max() * 1.05)
    ax.set_xlim(d.labtime.min(), d.labtime.max())

    text = "Records: " + "{:.2e}".format(col.attrs["Records"])
    bbox = dict(boxstyle="round", fc="blanchedalmond", ec="orange", alpha=0.5)
    ax.text(
        0.99,
        0.96,
        text,
        horizontalalignment="right",
        verticalalignment="top",
        transform=ax.transAxes,
        bbox=bbox,
    )


def PL_fig_plot(col, fitting=1):
    fig, gs = wt.artists.create_figure(
        width="double", nrows=3, default_aspect=0.25, hspace=0.7
    )
    axs = [plt.subplot(gs[i]) for i in range(3)]
    ylabels = [
        "$\\mathsf{counts \\; [Hz]}$",
        "$\\mathsf{counts}$",
        "$\\mathsf{cross-correlation}$",
    ]
    xlabels = [
        "$\\mathsf{labtime \\; [s]}$",
        "$\\mathsf{delay \\; [ns]}$",
        "$\\mathsf{\\tau/\\tau_{rep}}$",
    ]
    for ax, xlabel, ylabel in zip(axs, xlabels, ylabels):
        ax.grid()
        wt.artists.set_ax_labels(ax=ax, xlabel=xlabel, ylabel=ylabel, label_fontsize=22)
    PL_macrotime_plot(axs[0], col.macrohist, col)
    PL_picotime_plot(axs[1], col.picohist, fitting)
    PL_g2_plot(axs[2], col.g2hist)
    axs[0].set_title(col.attrs["identifier"])
    return fig, gs


def spectra_fig_plot(d):
    fig, gs = wt.artists.create_figure(default_aspect=0.5)
    ax = plt.subplot(gs[0])
    ax.plot(d, linewidth=2)
    wt.artists.set_ax_labels(xlabel=d.axes[0].label, ylabel="counts")
    ax.set_title(d.natural_name)
    ax.set_xlim(d.axes[0].min(), d.axes[0].max())
    ax.grid()
    return fig, gs


def confocal_scan_plot(d, sideplot=False):
    # d = d.copy()
    # d.create_channel(name='trace_log', values=np.sqrt(d['trace'][:]))
    # d.create_channel(name='retrace_log', values=np.sqrt(d['retrace'][:]))
    vmin = min(d["trace"].min(), d["retrace"].min())
    vmax = max(d["trace"].max(), d["retrace"].max())
    fig, gs = wt.artists.create_figure(width="double", cols=[1, 1, "cbar"])
    axs = [plt.subplot(gs[i]) for i in range(2)]
    axs[0].pcolor(d, channel="trace", vmin=vmin, vmax=vmax)
    axs[0].set_title("trace")
    axs[1].pcolor(d, channel="retrace", vmin=vmin, vmax=vmax)
    axs[1].set_title("retrace")
    wt.artists.set_fig_labels(
        xlabel="$\\mathsf{x \\; (\\mu m)}$",
        ylabel="$\\mathsf{y \\; (\\mu m)}$",
        title=d.natural_name,
    )
    cax = plt.subplot(gs[-1])
    ticks = np.linspace(vmin, vmax, 2)
    wt.artists.plot_colorbar(cax, ticks=ticks)
    if sideplot:
        first = [d.x.points, d.y.points, d.trace.points.T]
        second = [d.x.points, d.y.points, d.retrace.points.T]
        for ax, arrs_to_bin in zip(axs, [first, second]):
            wt.artists.add_sideplot(ax, along="x", arrs_to_bin=arrs_to_bin)
            wt.artists.add_sideplot(ax, along="y", arrs_to_bin=arrs_to_bin)

    return fig, gs
