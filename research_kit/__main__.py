import os

from .quick import directory_hl3_read_plot_save
from .quick import directory_princeton_spectrograph_ASCII_read_plot_save


def main():
    directory_hl3_read_plot_save(os.getcwd())
    directory_princeton_spectrograph_ASCII_read_plot_save(os.getcwd())


if __name__ == "__main__":
    main()
