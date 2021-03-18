import os

import argparse

from .quick import directory_hl3_read_plot_save
from .quick import directory_princeton_spectrograph_ASCII_read_plot_save
from .quick import directory_homebuilt_APDscan_ASCII_read_plot_save
from .quick import director_hl3_full


def read_plot_save():
    directory_hl3_read_plot_save(os.getcwd())
    directory_princeton_spectrograph_ASCII_read_plot_save(os.getcwd())
    directory_homebuilt_APDscan_ASCII_read_plot_save(os.getcwd())


def dir_hl3():
    parser = argparse.ArgumentParser(
        description="Crawl a directory and workup binary hl3 files to txt files"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Verbose? False by default"
    )
    parser.add_argument(
        "--hdf5",
        "-h5",
        action="store_true",
        help="store binary hd5 files False by default",
    )
    parser.add_argument(
        "--plot", "-p", action="store_true", help="plot results from each file"
    )
    parser.add_argument(
        "--fit", "-f", type=int, help="fit to a sum of exponentials: 1,2,3"
    )
    args = parser.parse_args()
    verbose = args.verbose
    hdf5 = args.hdf5
    plot = args.plot
    fit = args.fit
    directory = os.getcwd()
    director_hl3_full(directory, verbose, hdf5, plot, fit)
