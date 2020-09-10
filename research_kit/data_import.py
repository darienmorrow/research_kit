import numpy as np
import os
import struct
import pathlib
import warnings
import WrightTools as wt


def from_hl3(
    filepath,
    name=None,
    parent=None,
    keep_all=False,
    bin_which="zero",
    picobins=None,
    g2bins=None,
):
    filestr = os.fspath(filepath)
    filepath = pathlib.Path(filepath)
    if not ".hl3" in filepath.suffixes:
        warning.warn("wrong filetype. expected .hl3 file type")
    bytesize = os.path.getsize(filepath)
    ds = np.DataSource(None)
    f = ds.open(filestr, "rb")
    # unpack, decode, and calculate acquisition parameters
    acq_params = dict()
    places = [4, 33, 62, 91, 120]
    keys = ["TimeStamp", "DLLversion", "InitStat", "serial", "CalibStat"]
    for place, key in zip(places, keys):
        f.seek(place)
        acq_params[key] = f.read(25).strip().decode("utf-8")
    f.seek(145)
    keys = [
        "AcquisitionTime",
        "SyncRate",
        "CntRate0",
        "CntRate1",
        "Resolution",
        "SyncDivider",
        "SyncCFDLevel",
        "SyncZeroCross",
        "SyncOffset",
        "CFDLevel0",
        "CFDZeroCross0",
        "CFDOffset0",
        "CFDLevel1",
        "CFDZeroCross1",
        "CFDOffset1",
    ]
    for key in keys:
        acq_params[key] = struct.unpack("I", f.read(4))[0]
    acq_params["Records"] = (bytesize - f.tell()) // 4
    acq_params["syncperiod"] = 1e9 / acq_params["SyncRate"]
    # read in records
    syncperiod = acq_params["syncperiod"]
    resolution = acq_params["Resolution"]
    arr = np.fromfile(f, dtype=np.uint32)
    f.close()
    channeli = np.bitwise_and((arr // 33554432), 63)
    overflow_index = np.nonzero(channeli != 63)[0]
    picotime = (arr[overflow_index] // 1024) & 32767
    picotime = picotime * resolution * 1e-12
    out1 = np.bitwise_and(arr[overflow_index], 1023)
    out2 = 1024 * (overflow_index - np.arange(0, overflow_index.size, 1))
    macrotime = (out1 + out2) * syncperiod * 1e-9 + picotime
    # separate out the channels
    channels = channeli[overflow_index]
    # marker = np.nonzero((channels >= 2) & (channels <= 15))[0]
    picotime0 = picotime[channels == 0]
    picotime1 = picotime[channels == 1]
    macrotime0 = macrotime[channels == 0]
    macrotime1 = macrotime[channels == 1]

    if picobins is None:
        picobins = int(picotime.max() / (resolution * 1e-12) / 5)
    x = np.linspace(0, picotime.max(), picobins)
    hist0, bin_edges = np.histogram(picotime0, bins=picobins, range=(x[0], x[-1]))
    hist1, bin_edges = np.histogram(picotime1, bins=picobins, range=(x[0], x[-1]))
    histx = np.diff(bin_edges) / 2 + bin_edges[:-1]
    # offset all channel1 macrotime and picotimes
    offset = (
        histx[np.argmax(np.correlate(hist0, hist1, mode="same"))]
        - (picotime.max() - picotime.min()) / 2
    )
    macrotime[channels == 1] += offset
    picotime[channels == 1] += offset
    macrotime1 += offset
    picotime1 += offset
    # recalculate histogram
    hist1, bin_edges = np.histogram(picotime1, bins=picobins, range=(x[0], x[-1]))
    # do psuedo g2 calculation
    diffchan = np.diff(channels.astype(np.int32))
    diffmacro = np.diff(macrotime)
    g2 = diffchan * diffmacro
    g2 = g2[g2 != 0]
    g2 /= syncperiod * 1e-9
    # calculate area ratio (3 periods on either side of 0)
    a_1 = ((-3.5 < g2) & (g2 < -0.5)).sum()
    a1 = ((0.5 < g2) & (g2 < 3.5)).sum()
    a0 = ((-0.5 < g2) & (g2 < 0.5)).sum()
    arearatio = a0 / ((a_1 + a1) / 6)
    # calculate histograms
    macrobins = int(macrotime.max())
    if g2bins is None:
        g2bins = acq_params["Records"] // 4000
    g2y, bin_edges = np.histogram(g2, bins=g2bins, range=(-3.5, 3.5))
    g2x = np.diff(bin_edges) / 2 + bin_edges[:-1]
    if bin_which == "both":
        macroy, bin_edges = np.histogram(macrotime, bins=macrobins)
        macrox = np.diff(bin_edges) / 2 + bin_edges[:-1]
        picoy, bin_edges = np.histogram(picotime, bins=picobins)
        picox = np.diff(bin_edges) / 2 + bin_edges[:-1]
    elif bin_which == "zero":
        macroy, bin_edges = np.histogram(macrotime0, bins=macrobins)
        macrox = np.diff(bin_edges) / 2 + bin_edges[:-1]
        picoy, picox = hist0, histx
    elif bin_which == "one":
        macroy, bin_edges = np.histogram(macrotime1, bins=macrobins)
        macrox = np.diff(bin_edges) / 2 + bin_edges[:-1]
        picoy, picox = hist1, histx
    # now package it all into a Collection
    col = wt.Collection(name=name, parent=parent)
    col.attrs["identifier"] = filepath.stem + "_AcqTime-" + acq_params["TimeStamp"]
    col.attrs.update(acq_params)

    d = wt.Data(parent=col, name="macrohist", source=filestr)
    d.create_channel(name="counts", values=macroy)
    d.create_variable(name="labtime", units="s_t", values=macrox)
    d.transform("labtime")
    d = wt.Data(parent=col, name="picohist", source=filestr)
    d.create_channel(name="counts", values=picoy)
    d.create_variable(name="delay", units="ns", values=picox * 1e9)
    d.transform("delay")
    d = wt.Data(parent=col, name="g2hist", source=filestr)
    d.create_channel(name="counts", values=g2y)
    d.create_variable(name="delay", values=g2x)
    d.transform("delay")
    d.attrs["arearatio"] = arearatio
    if keep_all:
        c = wt.Collection(name="full_arrays", parent=col)
        # macrotime
        d = wt.Data(parent=c, name="macrotime", source=filestr)
        d.create_channel(name="macrotime", units="s_t", values=macrotime)
        # channels
        d = wt.Data(parent=c, name="channels", source=filestr)
        d.create_channel(name="channels", values=channels)
        # picotime
        d = wt.Data(parent=c, name="picotime", source=filestr)
        d.create_channel(name="picotime", units="ns", values=picotime * 1e9)
        # g2
        d = wt.Data(parent=c, name="g2", source=filestr)
        d.create_channel(name="g2", values=g2)
    return col
