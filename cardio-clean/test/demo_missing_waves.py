# coding: utf-8

import sys

from matplotlib import pyplot as plt

from cardio_clean.wavdetect import find_points
from cardio_clean.metadata import *

from cardio_clean.qrsdetect import qrs_detection
from demo_preprocessing import ecgread


def plot_with_markup(s, fs, metadata, chan):

    plt.cla()

    pt_keys = {"q_pos": "r","r_pos": "b", "s_pos": "b", "p_pos": "y", "t_pos":
    "g"}

    left_p = max(0, metadata[0]["qrs_center"]*fs - 100)
    right_p = min(s.shape[0]-1, metadata[-1]["qrs_center"]*fs + 100)
    plt.plot(s[int(left_p):int(right_p),chan], "b")

    for ncycle, qrs in enumerate(metadata):

        for k in pt_keys:
            point = qrs[k][chan]
            if point is not None:
                plt.scatter(point, s[point, chan], c=pt_keys[k])

        # подсветка P-зубца
        lb = qrs["p_start"][chan]
        rb = qrs["p_end"][chan]
        if all((lb, rb)):
            plt.plot(np.arange(lb, rb), s[lb:rb, chan], "y")

        # подсветка T-зубца
        lb = qrs["t_start"][chan]
        rb = qrs["t_end"][chan]
        if all((lb, rb)):
            plt.plot(np.arange(lb, rb), s[lb:rb, chan], "g")

        lb = qrs["st_start"][chan]
        rb = qrs["st_end"][chan]
        if all((lb, rb)):
            plt.plot(np.arange(lb, rb), s[lb:rb, chan], "r")

    plt.xlim((left_p, right_p))
    plt.grid(True)
    plt.show()


def show_waves(filename):

    sig, header = ecgread(filename)
    fs = header["fs"]
    if fs != 250:
        print("Warning! fs != 250 Hz, fs={}".format(fs))

    s = sig[:,0]

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
        s,
        fs=header["fs"]
    )

    qrsTypes = {}
    stcount=0

    suspicious = []

    for ncycle, qrs in enumerate(metadata):

        qrstype = qrs["qrsType"]
        qrsTypes[qrstype] = qrsTypes.get(qrstype, 0) + 1

        lb = qrs["st_start"][0]
        rb = qrs["st_end"][0]
        if all((lb, rb)):
            stcount += 1
        else:
            suspicious.append(ncycle)

    missing_hrt = [i for i,x in enumerate(metadata) if x["heartrate"] is None]
    print("Heartrate missing in beats\n{}".format(missing_hrt))

    print(qrsTypes)

    stdur = [qrs["st_duration"][0] for qrs in metadata if
             qrs["st_duration"][0] and qrs["st_duration"][0] is not None]

    print("{} cycles, {} ST segments, avg. {} ms of {} items".format(
        len(metadata),
        stcount,
        np.mean(stdur),
        len(stdur)
    ))

    for sus in suspicious:

        print("display cycle {}".format(sus))
        i1 = max(0, sus-1)
        i2 = min(len(metadata), sus+2)

        plot_with_markup(sig, fs, metadata[i1:i2], 0)

        print("Show next?")
        if sys.stdin.read(1).strip() == "n":
            break

if __name__ == "__main__":

    # Rh2022 = qr, noise
    # Rh2021 - Rs, extracyc
    # Rh2024 - p(q)RsT
    # Rh2025 = rs
    # filename = "/Users/arseniy/SERDECH/data/ROXMINE/Rh2024"
    filename = "TestFromDcm.ecg"

    if len(sys.argv) > 1:
        filename = sys.argv[0]

    show_waves(filename)


