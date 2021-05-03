import numpy as np
import os
import struct
import pathlib
import warnings
import WrightTools as wt


def _strip_end(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[: len(text) - len(suffix)]


def from_Cary(p):
    pass


def from_Horiba_2D(ps):
    pass

def from_RandS_trace(filepath, name=None, parent=None, cast_mag_to_dB=True):
    def caster(arr):
        if cast_mag_to_dB:
            return 20*np.log10(arr)
        else:
            return arr
    filestr = os.fspath(filepath)
    filepath = pathlib.Path(filepath)
    if not ".csv" in filepath.suffixes:
        warnings.warn("wrong filetype. expected .csv file type")
    # parse name
    if name is None:
        name = filepath.stem
    arr = np.genfromtxt(filepath, delimiter=';', skip_header=3, unpack=True)
    kwargs = {"name": name, "kind": "NetworkAnalizer", "source": filestr}
    d = wt.Data(parent=parent, **kwargs)
    d.create_variable('freq', units='GHz', values=arr[0]/1e9)
    chan_names = ['S11_re', 'S11_im', 'S12_re', 'S12_im', 'S21_re', 'S21_im', 'S22_re', 'S22_im']
    for i,c in enumerate(chan_names):
        d.create_channel(name=c, values=arr[i+1])
    d.transform('freq')
    chan_names = ['S11_mag', 'S12_mag', 'S21_mag', 'S22_mag']
    for i, c in enumerate(chan_names):
        a1 = arr[2*i+1]
        a2 = arr[2*i+2]
        mag = np.sqrt(a1**2 + a2**2)
        d.create_channel(name=c, values=caster(mag))
    return d

def from_homebuilt_APDscan_ASCII_triplet(partialfilepath, name=None, parent=None):
    ends = ["_A_set.txt", "_B_Mtr.txt", "_C_Mrt.txt"]
    for end in ends:
        partialfilepath = _strip_end(partialfilepath, end)
    filestrs = [os.fspath(partialfilepath + end) for end in ends]
    fs = [np.DataSource(None).open(filestr, "rb") for filestr in filestrs]
    # first harvest some metadata from A_set file
    arr = np.loadtxt(fs[0], max_rows=1, dtype=int)
    x0, y0, extent, pixels = arr[0], arr[1], arr[2], arr[3]
    x = np.linspace(x0, extent + x0, pixels)
    y = np.linspace(y0, extent + y0, pixels)
    # grab trace and retrace data
    trace = np.genfromtxt(fs[1], unpack=True)
    retrace = np.genfromtxt(fs[2], unpack=True)
    # parse name
    if name is None:
        name = pathlib.Path(partialfilepath).stem
    # create data
    kwargs = {"name": name, "kind": "APDscan", "source": filestrs[1]}
    if parent is None:
        data = wt.Data(**kwargs)
    else:
        data = parent.create_data(**kwargs)
    data.create_variable("x", values=x[:, None], units="um")
    data.create_variable("y", values=y[None, :], units="um")
    data.create_channel("trace", values=trace)
    data.create_channel("retrace", values=retrace)
    data.transform("x", "y")
    for f in fs:
        f.close()
    return data


def from_princeton_spectrograph_ASCII(filepath, name=None, parent=None, keep_2D=False):
    filestr = os.fspath(filepath)
    filepath = pathlib.Path(filepath)
    if not ".txt" in filepath.suffixes:
        warnings.warn("wrong filetype. expected .txt file type")
    ds = np.DataSource(None)
    f = ds.open(filestr, "rb")
    # array creation
    # first we figure out what type of data format we are dealing with
    arr = np.genfromtxt(f, max_rows=1, unpack=True)
    f.seek(0)
    numnan = np.isnan(arr).sum()
    if numnan == 3:
        z = np.genfromtxt(f, dtype=int, skip_header=3, unpack=True)
        cols = tuple(np.arange(3, z.shape[0] + 1, 1, dtype=int))
        f.seek(0)
        wl = np.genfromtxt(f, usecols=cols, skip_header=1, skip_footer=z.shape[1] + 1)
        frame = z[0, :]
        strip = z[1, :]
        z = z[2:, :]
    elif numnan == 1:
        z = np.genfromtxt(f, dtype=int, skip_header=3, unpack=True)
        cols = tuple(np.arange(1, z.shape[0] + 1, 1, dtype=int))
        f.seek(0)
        wl = np.genfromtxt(f, skip_header=1, max_rows=1)
        strip = z[0, :]
        wl = wl[1:]
        z = z[1:, :]
    elif arr.size == 3:
        arrs = []
        dtypes = [float, int, int]
        for i, col in enumerate(range(3)):
            f.seek(0)
            arrs.append(np.genfromtxt(f, unpack=True, usecols=(col), dtype=dtypes[i]))
        length1 = arrs[1].max()
        length0 = arrs[1].shape[0] // length1
        newshape = (length1, length0)
        for i in range(len(arrs)):
            arrs[i] = np.reshape(arrs[i], newshape)
        wl = arrs[0][0, :]
        z = arrs[2].T
        strip = np.linspace(1, z.shape[-1], z.shape[-1], dtype=int)
    else:
        raise Exception("data format not currently supported")

    # parse name
    if name is None:
        name = filepath.stem
    # create data
    kwargs = {"name": name, "kind": "Princeton", "source": filestr}
    if parent is None:
        data = wt.Data(**kwargs)
    else:
        data = parent.create_data(**kwargs)
    if keep_2D:
        data.create_variable("wavelength", units="nm", values=wl[:, None])
        data.create_variable("strip", values=strip[None, :])
        data.create_channel("counts", values=z)
        data.transform("wavelength", "strip")
    else:
        z = np.sum(z, axis=1)
        data.create_variable("wavelength", units="nm", values=wl)
        data.create_channel("counts", values=z)
        data.transform("wavelength")
    f.close()
    return data


def from_hl3(
    filepath,
    name=None,
    parent=None,
    keep_all=False,
    bin_which="zero",
    picobins=None,
    g2bins=None,
    area_ratio_calc = True,
):
    filestr = os.fspath(filepath)
    filepath = pathlib.Path(filepath)
    if not ".hl3" in filepath.suffixes:
        warnings.warn("wrong filetype. expected .hl3 file type")
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
    macrotime = (out1 + out2) / acq_params["SyncRate"] + picotime
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
    if ((hist0.max() > 1) & (hist1.max() > 1)):
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
    _a = np.array([a_1, a1, a0])
    if np.any(_a == 0):
        arearatio = 0
    else:
        if area_ratio_calc:
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
    if bin_which == "combined":
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
    elif bin_which == "both":
        macroy0, bin_edges = np.histogram(macrotime0, bins=macrobins)
        macroy1, bin_edges = np.histogram(macrotime1, bins=macrobins)
        macrox = np.diff(bin_edges) / 2 + bin_edges[:-1]
        picoy0, picoy1, picox = hist0, hist1, histx
    # now package it all into a Collection
    if name is None:
        name = filepath.stem
    col = wt.Collection(name=name, parent=parent)
    col.attrs["identifier"] = filepath.stem + "_AcqTime-" + acq_params["TimeStamp"]
    col.attrs.update(acq_params)
    if bin_which == "both":
        d = wt.Data(parent=col, name="macrohist", source=filestr)
        d.create_channel(name="counts0", values=macroy0)
        d.create_channel(name="counts1", values=macroy1)
        d.create_variable(name="labtime", units="s_t", values=macrox)
        d.transform("labtime")
        d = wt.Data(parent=col, name="picohist", source=filestr)
        d.create_channel(name="counts0", values=picoy0)
        d.create_channel(name="counts1", values=picoy1)
        d.create_variable(name="delay", units="ns", values=picox * 1e9)
        d.transform("delay")
        d = wt.Data(parent=col, name="g2hist", source=filestr)
        d.create_channel(name="counts", values=g2y)
        d.create_variable(name="delay", values=g2x)
        d.transform("delay")
    else:
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
    if area_ratio_calc:
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
