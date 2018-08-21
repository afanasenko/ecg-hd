# coding: utf-8

import sys

from matplotlib import pyplot as plt

from cardio_clean.wavdetect import find_points
from cardio_clean.metadata import *

from cardio_clean.sigbind import fix_baseline
from cardio_clean.qrsdetect import qrs_detection
from cardio_clean.util import ecgread


def plot_with_markup(sig, fs, metadata):

    numch = len(metadata[0]["r_pos"])

    pt_keys = {"q_pos": "b","r_pos": "r", "s_pos": "b", "p_pos": "y", "t_pos":
    "g"}

    left_p = max(0, metadata[0]["qrs_center"]*fs - 120)
    right_p = min(sig.shape[0]-1, metadata[-1]["qrs_center"]*fs + 120)

    tscale = 1.0 / fs

    f, ax = plt.subplots(numch, 1, sharex="col")

    for chan, s in signal_channels(sig):
        lb = int(left_p)
        rb = int(right_p)
        ax[chan].plot(np.arange(lb, rb)*tscale, s[lb:rb], "b")

        for ncycle, qrs in enumerate(metadata):

            for k in pt_keys:
                point = qrs[k][chan]
                if point is not None:
                    ax[chan].scatter(point*tscale, s[point], c=pt_keys[k])

            # подсветка R-зубца
            lb = qrs["r_start"][chan]
            rb = qrs["r_end"][chan]
            if all((lb, rb)):
                ax[chan].plot(np.arange(lb, rb) * tscale, s[lb:rb], "r")

            # подсветка P-зубца
            lb = qrs["p_start"][chan]
            rb = qrs["p_end"][chan]
            if all((lb, rb)):
                ax[chan].plot(np.arange(lb, rb)*tscale, s[lb:rb], "y")

            # подсветка T-зубца
            lb = qrs["t_start"][chan]
            rb = qrs["t_end"][chan]
            if all((lb, rb)):
                ax[chan].plot(np.arange(lb, rb)*tscale, s[lb:rb], "g")

            lb = qrs["st_start"][chan]
            rb = qrs["st_end"][chan]
            if all((lb, rb)):
                ax[chan].plot(np.arange(lb, rb)*tscale, s[lb:rb], "m")

    plt.show()


def show_waves(filename):

    sig, header = ecgread(filename)
    fs = header["fs"]
    if fs != 250:
        print("Warning! fs != 250 Hz, fs={}".format(fs))

    sig = fix_baseline(
        sig,
        fs=fs,
        bias_window_ms=1500
    )

    metadata = qrs_detection(
        sig[:,:],
        fs=header["fs"]
    )[0]

    find_points(
        sig[:, :],
        fs=header["fs"],
        bias=header["baseline"],
        metadata=metadata
    )

    metadata_postprocessing(
        metadata,
        sig,
        header
    )

    suspicious = []

    for ncycle, qrs in enumerate(metadata):

        if max(qrs["qt_duration"]) > 1.0:
            suspicious.append(ncycle)

    missing_hrt = [i for i,x in enumerate(metadata) if x["heartrate"] is None]
    print("Heartrate missing in beats\n{}".format(missing_hrt))

    for sus in suspicious:

        print("display cycle {}, QT = {}".format(
            sus,
            metadata[sus]["qt_duration"]
        ))
        i1 = sus
        i2 = min(len(metadata), sus+2)

        plot_with_markup(sig, fs, metadata[i1:i2])
        # plot_with_markup(sig, fs, [metadata[sus]])

        print("Show next?")
        if sys.stdin.read(1).strip() == "n":
            break

if __name__ == "__main__":

    # Rh2022 = qr, noise
    # Rh2021 - Rs, extracyc
    # Rh2024 - p(q)RsT
    # Rh2025 = rs
    filename = "/Users/arseniy/SERDECH/data/ROXMINE/Rh2024"
    # filename = "TestFromDcm.ecg"

    if len(sys.argv) > 1:
        filename = sys.argv[0]

    show_waves(filename)


