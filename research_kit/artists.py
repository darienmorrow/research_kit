import WrightTools as wt
import numpy as np
import matplotlib.pyplot as plt

from . import fit as fit


def PL_g2_plot(ax, d):
    ax.plot(d, linewidth=1)
    ax.fill_between(d.delay.points, d.counts.points, alpha=0.3)
    ax.set_xlim(d.delay.min(), d.delay.max())
    ax.set_ylim(0, d.counts.max() * 1.05)
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


def PL_picotime_plot(ax, d, fitting=True):
    x = d.delay.points
    y = d.counts.points
    ax.plot(d, linewidth=1, alpha=1)
    ax.fill_between(x, y, alpha=0.3)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min() + 1, y.max() * 1.1)
    ax.set_yscale("log")
    if fitting:
        pfit, perr, ymodel = fit.exp_fit(x, y)
        ax.plot(x, ymodel, color="C1", linewidth=4, alpha=0.75)
        t1 = pfit[0]
        text = "$\\mathsf{\\tau_1 =" + str(round(np.abs(t1), 1)) + "\\;ns}$"
        bbox = dict(boxstyle="round", fc="C1", ec="C1", alpha=0.5)
        ax.text(
            0.99,
            0.96,
            text,
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
            bbox=bbox,
        )
        pfit, perr, ymodel = fit.biexp_fit(x, y)
        ax.plot(x, ymodel, color="C2", linewidth=4, alpha=0.75)
        text = (
            "$\\mathsf{\\tau_1 ="
            + str(round(np.abs(pfit[0]), 1))
            + "\\;ns, \\; \\tau_2 ="
            + str(round(np.abs(pfit[2]), 1))
            + "\\;ns, \\;"
            + "A_2/A_1 ="
            + str(round(np.abs(pfit[3] / pfit[1]), 2))
            + "}$"
        )
        bbox = dict(boxstyle="round", fc="C2", ec="C12", alpha=0.5)
        ax.text(
            0.99,
            0.8,
            text,
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
            bbox=bbox,
        )


def PL_macrotime_plot(ax, d, col):
    ax.plot(d, linewidth=1)
    ax.fill_between(d.labtime.points, d.counts.points, alpha=0.3)
    ax.set_xlim(d.labtime.min(), d.labtime.max())
    ax.set_ylim(d.counts.min() + 1, d.counts.max() * 1.05)
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


def PL_fig_plot(col, fitting=True):
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
        wt.artists.set_ax_labels(ax=ax, xlabel=xlabel, ylabel=ylabel)
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
