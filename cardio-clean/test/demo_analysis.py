# coding: utf-8

import time
from matplotlib import pyplot as plt
from scipy.signal import hann
from scipy.fftpack import fft

from cardio_clean.wavdetect import ddwt, find_points, zcfind
from cardio_clean.arrythmia import *
from cardio_clean.ishemia import mock_ishemia_episodes, define_ishemia_episodes

from cardio_clean.sigbind import fix_baseline
from cardio_clean.qrsdetect import qrs_detection
from cardio_clean.util import ecgread

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


def show_decomposition(filename, ch, lim):

    sig, hdr = ecgread(filename)

    if lim:
        lim = min(lim, sig.shape[0])
    else:
        lim = sig.shape[0]
    s = sig[:lim, ch]

    app, ders = ddwt(s, num_scales=5)
    multiplot(ders)


def print_summary(metadata, ch=1):
    classes = {}
    qrs_types = {}
    st_intervals = []
    rr_intervals = []
    wcount = {w: 0 for w in ("p", "q", "r", "s", "t")}

    for ncycle, qrs in enumerate(metadata):
        qrstype = qrs["qrsType"]
        qrs_types[qrstype] = qrs_types.get(qrstype, 0) + 1

        classletter = qrs["complex_type"]
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

        if qrs["st_duration"][ch] is not None:
            st_intervals.append(qrs["st_duration"][ch])

        if qrs["RR"] is not None:
            rr_intervals.append(qrs["RR"])

    print("Комплексы: {}".format(len(metadata)))
    print("Зубцы: {}".format(wcount))

    print("Конфигурации qrs:")
    print(qrs_types)
    print("Типы комплексов:")
    print(classes)

    print("ST-интервалы: {}, в среднем {} мс".format(
        len(st_intervals), np.mean(st_intervals)
    ))

    print("RR-интервалы: {}, в среднем {} мс".format(
        len(rr_intervals), np.mean(rr_intervals)
    ))


def show_qrs(filename, chan, lim):
    sig, header = ecgread(filename)

    fs = header["fs"]
    if fs != 250:
        print("Warning! fs={}".format(fs))

    if lim:
        lim = min(lim, sig.shape[0])
    else:
        lim = sig.shape[0]

    metadata, pant = qrs_detection(
        sig[:lim,:],
        fs=header["fs"]
    )

    approx, ders = ddwt(sig[:lim, chan], num_scales=5)

    fig, axarr = plt.subplots(len(ders), 1, sharex="col")
    for i, s in enumerate(ders):
        axarr[i].plot(s,"b")
        axarr[i].plot(approx[i],"g")

        for ncycle, qrs in enumerate(metadata):
            lb = int(qrs["qrs_start"] * fs)
            rb = int(qrs["qrs_end"] * fs)
            axarr[i].plot(np.arange(lb, rb), s[lb:rb], "r")

        axarr[i].grid()

    ms_step = 1000
    locs = np.arange(0, len(ders[0]), ms_to_samples(ms_step, fs))
    labs = [time.strftime("%M:%S", time.gmtime(x * ms_step)) for x in locs]
    axarr[-1].set_xticks(locs, labs)
    ms_step = 100
    locs = np.arange(0, len(ders[0]), ms_to_samples(ms_step, fs))
    axarr[-1].set_xticks(locs, minor=True)

    print("Look at the plots")
    plt.show()


def show_waves(filename, chan, lim):
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
        metadata=metadata
    )

    metadata_postprocessing(
        metadata,
        sig[:lim, :],
        header
    )

    plt.plot(s, "b")

    # Цвета для раскрашивания зубцов на графике
    pt_keys = {"q_pos": "r","r_pos": "g", "s_pos": "b", "p_pos": "y", "t_pos":
    "m"}

    for ncycle, qrs in enumerate(metadata):

        lb = int(qrs["qrs_start"] * fs)
        rb = int(qrs["qrs_end"] * fs)
        iso = qrs["isolevel"][chan] + header["baseline"][chan]
        plt.plot([lb, rb], [iso]*2, "g:")

        for k in pt_keys:
            point = qrs[k][chan]
            if point is not None:
                plt.scatter(point, s[point], c=pt_keys[k])

        if qrs["qrs_center"] > 60:
            break

        lb = qrs["r_start"][chan]
        rb = qrs["r_end"][chan]
        if all((lb, rb)):
            plt.plot(np.arange(lb, rb), s[lb:rb], "g")

        lb = qrs["st_start"][chan]
        rb = qrs["st_end"][chan]
        if all((lb, rb)):
            plt.plot(np.arange(lb, rb), s[lb:rb], "r")

        if qrs["complex_type"] == "V":

            if qrs["r_pos"][chan] is not None:
                plt.text(
                    qrs["r_pos"][chan],
                    qrs["r_height"][chan],
                    qrs["complex_type"],
                    color="r"
                )

        plt.text(
            qrs["qrs_start"]*fs,
            np.min(s),
            str(ncycle)
        )

    #missing_hrt = [i for i,x in enumerate(metadata) if x["heartrate"] is None]
    #print("Heartrate missing in beats\n{}".format(missing_hrt))

    #plt.xlim((200,700))

    #r = define_rythm(metadata)
    #print(r)

    print("Ишемия...")
    m = define_ishemia_episodes(sig[:lim, :], header, metadata)
    print(m)
    print("Число эпизодов: {}".format(len(m)))

    #for x in r:
    #    print("ритм {}: {} с".format(x["desc"], x["end"] - x["start"]))

    print_summary(metadata)

    plt.grid(True)
    plt.show()

if __name__ == "__main__":

    # Rh2022 = qr, noise
    # Rh2021 - Rs, extracyc
    # Rh2024 - p(q)RsT
    # Rh2025 = rs
    # Rh2010 - дрейф, шум, артефакты

    #filename = "/Users/arseniy/SERDECH/data/PHYSIONET/I60"
    filename = "TestFromDcm.ecg"

    #show_filter_responses()
    #show_decomposition(filename, 1, 50000)
    #show_qrs(filename, 1, 20000)
    show_waves(filename, 1, 20000)

