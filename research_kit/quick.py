import os
import WrightTools as wt
import matplotlib.pyplot as plt
from itertools import islice
import numpy as np

from . import artists as art
from . import data_import as di


def single_homebuilt_APDscan_read_plot_save(
    partialfilepath, overwrite=True, verbose=False, sideplot=False
):
    d = di.from_homebuilt_APDscan_ASCII_triplet(partialfilepath)
    fig, gs = art.confocal_scan_plot(d, sideplot=sideplot)
    p = partialfilepath + ".wt5"
    d.save(p, overwrite=overwrite, verbose=verbose)
    p = partialfilepath + ".png"
    if (os.path.exists(p) == False) or overwrite:
        wt.artists.savefig(p, fig=fig)


def directory_homebuilt_APDscan_ASCII_read_plot_save(
    directory, overwrite=True, verbose=False
):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith("_A_set.txt"):
                p = os.path.join(root, file)
                # Do a crude test to make sure the txt file isn't a regular file.
                go = False
                with open(p) as lines:
                    for line in islice(lines, 4, 5):
                        if line == "1APD":
                            go = True
                if go:
                    single_homebuilt_APDscan_read_plot_save(
                        p, overwrite=overwrite, verbose=verbose
                    )


def single_hl3_read_plot_save(filepath, overwrite=True, verbose=False, fitting=True):
    col = di.from_hl3(filepath)
    fig, gs = art.PL_fig_plot(col, fitting=fitting)

    pre_name = os.path.splitext(filepath)[0]
    p = pre_name + ".wt5"
    col.save(p, overwrite=overwrite, verbose=verbose)
    p = pre_name + ".png"
    if (os.path.exists(p) == False) or overwrite:
        wt.artists.savefig(p, fig=fig)

def directory_hl3_ASCII_save(directory, overwrite=True, verbose=False):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".hl3"):
                p = os.path.join(root, file)
                col = di.from_hl3(p, bin_which = "both")
                d = col.picohist
                pre_name = os.path.splitext(p)[0]
                outp = pre_name + '_hl3_export.txt'
                x = d.delay.points
                y0 = d.counts0.points
                y1 = d.counts1.points
                arr = np.stack((x, y0, y1), axis=-1)
                if (os.path.exists(outp) == False) or overwrite:
                    np.savetxt(outp, arr, delimiter=",")
                

def directory_hl3_read_plot_save(
    directory, overwrite=True, verbose=False, fitting=True
):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".hl3"):
                p = os.path.join(root, file)
                single_hl3_read_plot_save(
                    p, overwrite=overwrite, verbose=verbose, fitting=fitting
                )


def single_princeton_spectrograph_ASCII_read_plot_save(
    filepath, overwrite=True, verbose=False
):
    d = di.from_princeton_spectrograph_ASCII(filepath)
    fig, gs = art.spectra_fig_plot(d)
    pre_name = os.path.splitext(filepath)[0]
    p = pre_name + ".wt5"
    d.save(p, overwrite=overwrite, verbose=verbose)
    p = pre_name + ".png"
    if (os.path.exists(p) == False) or overwrite:
        plt.savefig(p, dpi=300, bbox_inches="tight")
        plt.close()


def directory_princeton_spectrograph_ASCII_read_plot_save(
    directory, overwrite=True, verbose=False
):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".txt"):
                p = os.path.join(root, file)
                # do a crude test to make sure the txt file isn't a regular file.
                go = False
                with open(p) as lines:
                    for line in islice(lines, 2, 3):
                        if line[1:6] == "Frame":
                            go = True
                if go:
                    single_princeton_spectrograph_ASCII_read_plot_save(
                        p, overwrite=overwrite, verbose=verbose
                    )
