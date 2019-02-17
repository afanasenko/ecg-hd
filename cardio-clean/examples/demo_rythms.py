# coding: utf-8

import time
from matplotlib import pyplot as plt

from cardio_clean.wavdetect import find_points
from cardio_clean.arrythmia import *

from cardio_clean.sigbind import fix_baseline
from cardio_clean.qrsdetect import qrs_detection
from cardio_clean.util import ecgread


def scatterogram(series, ax):

    for i, x in enumerate(series):
        if i:
            ax.scatter(series[i-1], x, c="k")


def show_qt(filename, chan, lim):
    sig, header = ecgread(filename)

    fs = header["fs"]
    if fs != 250:
        print("Warning! fs={}".format(fs))

    print("Усиление: {}".format(header["adc_gain"]))

    sig = fix_baseline(
        sig,
        fs=fs,
        bias_window_ms=1500
    )

    if lim:
        lim = min(lim, sig.shape[0])
    else:
        lim = sig.shape[0]
    s = sig[:lim, chan]

    metadata, foo = qrs_detection(
        sig[:lim,:],
        fs=header["fs"]
    )

    find_points(
        sig[:lim, :],
        fs=header["fs"],
        bias=header["baseline"],
        gain=header["adc_gain"],
        metadata=metadata
    )

    metadata_postprocessing(
        metadata,
        sig[:lim, :],
        header
    )

    qt_int = []
    qtc_int = []
    rr_int = []


    for ncycle, qrs in enumerate(metadata):

        #if is_artifact(qrs):
        #    continue

        if qrs["qt_duration"] is not None:
            qt_int.append(qrs["qt_duration"])

        if qrs["qtc_duration"] is not None:
            qtc_int.append(qrs["qtc_duration"])

        if qrs["RR"] is not None:
            rr_int.append(qrs["RR"])

    fig, axarr = plt.subplots(2, 2)

    axarr[0,0].set_title("QT")
    scatterogram(qt_int, axarr[0,0])
    axarr[0, 1].set_title("QTc")
    scatterogram(qtc_int, axarr[0,1])
    axarr[1, 0].set_title("RR")
    scatterogram(rr_int, axarr[1,0])

    plt.show()

if __name__ == "__main__":

    filename = "/Users/arseniy/SERDECH/data/PHYSIONET/I60"

    show_qt(filename, 1, 20000)

