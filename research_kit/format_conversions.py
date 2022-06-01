import numpy as np
import tidy_headers
import collections

def arr_to_ASCII(arr, column_names, pout):
    """
    Write array to ASCII file.

    Parameters
    ----------
    arr : array_like
        vstack is useful here arr = np.vstack(outs).T
    column_names : list
        names of columns
    outp : str
        filepath to write files to

    Returns
    -------
    None.

    """
    # vstack is useful here arr = np.vstack(outs).T
    out_dict = collections.OrderedDict()
    out_dict["columns"] = column_names
    tidy_headers.write(pout, out_dict)

    with open(pout, "ab") as f:
        np.savetxt(f, arr, delimiter="\t")


def col_TRPL_to_ASCII(col, outp):
    """
    Convert a wt5 TRPL collection to 3 ASCII files.

    Parameters
    ----------
    col : wt5 collection
        collection of TRPL data from hl3 file
    outp : str
        filepath to write files to. Will have names and .txt appended

    Returns
    -------
    None.

    """
    metakeys = [
        "TimeStamp",
        "Records",
        "Resolution",
        "syncperiod",
        "AcquisitionTime",
        "name",
    ]
    metavals = [col.attrs[key] for key in metakeys]
    meta = collections.OrderedDict()
    for key, val in zip(metakeys, metavals):
        meta[key] = val
    for name in col:
        d = col[name]
        pout = str(outp) + "_" + name + ".txt"
        outs = []
        outs_names = []
        for axe in d.axes:
            outs.append(axe.points)
            if axe.units == None:
                n = axe.natural_name
            else:
                n = axe.natural_name + " [" + axe.units + "]"
            outs_names.append(n)
        for chan in d.channels:
            outs.append(chan.points)
            outs_names.append(chan.natural_name)
        out_arr = np.vstack(outs).T
        out_dict = meta.copy()

        out_dict["columns"] = outs_names
        tidy_headers.write(pout, out_dict)

        with open(pout, "ab") as f:
            np.savetxt(f, out_arr, delimiter="\t")
