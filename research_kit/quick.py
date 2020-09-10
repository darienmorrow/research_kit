import os
import WrightTools as wt

from . import artists as art
from . import data_import as di



def single_hl3_read_plot_save(filepath, overwrite=True, verbose=False, fitting=True):
    col = di.from_hl3(filepath)
    fig, gs = art.PL_fig_plot(col, fitting=fitting)
    
    pre_name = os.path.splitext(filepath)[0]
    p = pre_name + '.wt5'
    col.save(p, overwrite=overwrite, verbose=verbose)
    p = pre_name + '.png'
    if (os.path.exists(p) == False) or overwrite:
        wt.artists.savefig(p, fig=fig)
        
def directory_hl3_read_plot_save(directory, overwrite=True, verbose=False, fitting=True):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".hl3"):
                 p = os.path.join(root, file)
                 single_hl3_read_plot_save(p, overwrite=overwrite, verbose=verbose, fitting=fitting)
        

