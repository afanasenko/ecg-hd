# coding: utf-8

from matplotlib import pyplot as plt
from scipy.signal import lfilter, hann
from scipy.fftpack import fft

from cardio_clean.wavdetect import ddwt, find_points, zcfind
from cardio_clean.metadata import *

from cardio_clean.qrsdetect import qrs_detection
from demo_preprocessing import ecgread


def dummy_shift(x, n):
    return np.concatenate((x[n:], np.ones(n)*x[-1]))


def multiplot(siglist):
    fig, axarr = plt.subplots(len(siglist), 1, sharex="col")
    for i, s in enumerate(siglist):
        axarr[i].plot(s)
        axarr[i].grid()
    print("Look at the plots")
    plt.show()


def multispectrum(siglist):
    fig, axarr = plt.subplots(len(siglist), 1, sharex="col")
    for i, s in enumerate(siglist):
        n1 = len(s)
        n2 = int(n1/2)
        wgt = np.array(hann(n1))
        yf = fft(s * wgt)
        yf = np.abs(yf[0:n2])
        axarr[i].plot(yf)
        axarr[i].grid()
    print("Look at the plots")
    plt.show()


def multilinespec(siglist):
    leg = []
    for i, s in enumerate(siglist):
        n1 = len(s)
        n2 = int(n1/2)
        wgt = np.array(hann(n1))
        yf = fft(s * wgt)
        yf = np.abs(yf[0:n2])
        plt.plot(yf)
        leg.append(str(i))
    plt.grid()
    plt.legend(leg)
    print("Look at the plots")
    plt.show()


def show_filter_responses(spectral=False):

    T = 256
    x = np.zeros(T, float)
    mid = int(T/2)
    x[mid] = 1.0
    app, ders = ddwt(x, num_scales=5)

    if spectral:
        multilinespec(ders[1:])
    else:
        multiplot(ders)
        for i, x in enumerate(ders[1:]):
            dt_theor = (2.0**i - 1)/2
            dt_exrep = zcfind(x) - mid
            print(dt_theor, dt_exrep)


def show_decomposition():

    sig, hdr = ecgread("/Users/arseniy/SERDECH/data/ROXMINE/Rh2024")
    app, ders = ddwt(sig[:20000,0], num_scales=5)
    multiplot(ders)


def print_summary(metadata, ch=1):
    classes = {}
    qrs_types = {}

    wcount = {w: 0 for w in ("p", "q", "r", "s", "t")}

    for ncycle, qrs in enumerate(metadata):
        qrstype = qrs["qrsType"]
        qrs_types[qrstype] = qrs_types.get(qrstype, 0) + 1

        classletter = qrs["qrs_class_id"][0]
        classes[classletter] = classes.get(classletter, 0) + 1

        if qrs["p_pos"][ch] is not None:
            wcount["p"] += 1
        if qrs["q_pos"][ch] is not None:
            wcount["q"] += 1
        if qrs["r_pos"][ch] is not None:
            wcount["r"] += 1
        if qrs["s_pos"][ch] is not None:
            wcount["s"] += 1
        if qrs["t_pos"][ch] is not None:
            wcount["t"] += 1

    print("Зубцы:")
    print(wcount)
    print("Конфигурации qrs:")
    print(qrs_types)
    print("Типы комплексов:")
    print(classes)


def show_waves():
    # Rh2022 = qr, noise
    # Rh2021 - Rs, extracyc
    # Rh2024 - p(q)RsT
    # Rh2025 = rs
    sig, header = ecgread("/Users/arseniy/SERDECH/data/ROXMINE/Rh2025")
    #sig, header = ecgread("TestFromDcm.ecg")
    fs = header["fs"]
    if fs != 250:
        print("Warning! fs={}".format(fs))

    chan = 1
    lim = min(20000000, sig.shape[0])
    s = sig[:lim,chan]

    metadata, foo = qrs_detection(
        sig[:lim,:],
        fs=header["fs"]
    )

    find_points(
        sig[:lim, :],
        fs=header["fs"],
        bias=header["baseline"],
        metadata=metadata
    )

    metadata_postprocessing(
        metadata,
        sig[:lim, :],
        fs=header["fs"]
    )

    plt.plot(s, "b")

    stcount=0

    pt_keys = {"q_pos": "r","r_pos": "b", "s_pos": "b", "p_pos": "y", "t_pos":
    "g"}

    for ncycle, qrs in enumerate(metadata):

        lb = int(qrs["qrs_start"] * fs)
        rb = int(qrs["qrs_end"] * fs)
        iso = qrs["isolevel"][chan] + header["baseline"][chan]
        plt.plot([lb, rb], [iso]*2, "g:")

        for k in pt_keys:
            point = qrs[k][chan]
            if point is not None and point < 15500:
                plt.scatter(point, s[point], c=pt_keys[k])

        lb = qrs["st_start"][chan]
        rb = qrs["st_end"][chan]
        if all((lb, rb)):
            plt.plot(np.arange(lb, rb), s[lb:rb], "r")
            stcount += 1

    missing_hrt = [i for i,x in enumerate(metadata) if x["heartrate"] is None]
    print("Heartrate missing in beats\n{}".format(missing_hrt))

    plt.xlim((200,700))

    stdur = [qrs["st_duration"][chan] for qrs in metadata if
             qrs["st_duration"][chan]]

    print("{} cycles, {} ST segments, avg. {} ms of {} cycles".format(
        len(metadata),
        stcount,
        np.mean(stdur) if stdur else "-",
        len(stdur)
    ))

    print_summary(metadata)
    plt.show()

if __name__ == "__main__":

    #show_filter_responses()
    #show_decomposition()
    show_waves()

