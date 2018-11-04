# coding: utf-8

import time
import os
from matplotlib import pyplot as plt
from scipy.signal import hann
from scipy.fftpack import fft

from cardio_clean.wavdetect import ddwt, find_points, zcfind
from cardio_clean.arrythmia import *
from cardio_clean.ishemia import mock_ishemia_episodes, define_ishemia_episodes

from cardio_clean.sigbind import fix_baseline
from cardio_clean.qrsdetect import qrs_detection
from cardio_clean.util import ecgread, common_signal_names


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


def show_decomposition(filename, chan, smp_from=0, smp_to=0):

    sig, hdr = ecgread(filename)

    if smp_to:
        smp_to = min(smp_to, sig.shape[0])
    else:
        smp_to = sig.shape[0]
    s = sig[smp_from:smp_to, chan]

    app, ders = ddwt(s, num_scales=5)
    multiplot(ders)


def show_decomposition2(filename, chan, smp_from=0, smp_to=0):

    sig, hdr = ecgread(filename)

    if smp_to:
        smp_to = min(smp_to, sig.shape[0])
    else:
        smp_to = sig.shape[0]

    s = sig[smp_from:smp_to, chan]

    app, ders = ddwt(s, num_scales=5)
    # multiplot(ders)
    fig, axarr = plt.subplots(5, 1, sharex="col")
    for i in range(5):
        axarr[i].plot(s)
        if i>0:
            axarr[i].plot(app[i])
        axarr[i].grid()
    print("Look at the plots")
    plt.show()



def show_raw(filename, smp_from=0, smp_to=0):

    sig, header = ecgread(filename)
    fs = header["fs"]

    sig = fix_baseline(
        sig,
        fs=fs,
        bias_window_ms=1500
    )

    if smp_to:
        smp_to = min(smp_to, sig.shape[0])
    else:
        smp_to = sig.shape[0]

    fig, axarr = plt.subplots(sig.shape[1], 1, sharex="col")
    for i, s in signal_channels(sig):
        axarr[i].plot(s[smp_from:smp_to])
        axarr[i].grid()
    print("Look at the plots")
    plt.show()


def print_summary(metadata, chan):
    classes = {}
    classids = {}
    st_intervals = []
    qt_intervals = []
    rr_intervals = []
    wcount = {w: 0 for w in ("p", "q", "r", "s", "t")}

    for ncycle, qrs in enumerate(metadata):

        classletter = qrs["complex_type"]
        classes[classletter] = classes.get(classletter, 0) + 1

        classid = qrs["qrs_class_id"]
        classids[classid] = classids.get(classid, 0) + 1

        if qrs["p_pos"][chan] is not None:
            wcount["p"] += 1
        if qrs["q_pos"][chan] is not None:
            wcount["q"] += 1
        if qrs["r_pos"][chan] is not None:
            wcount["r"] += 1
        if qrs["s_pos"][chan] is not None:
            wcount["s"] += 1
        if qrs["t_pos"][chan] is not None:
            wcount["t"] += 1

        if qrs["st_duration"][chan] is not None:
            st_intervals.append(qrs["st_duration"][chan])

        if qrs["RR"] is not None:
            rr_intervals.append(qrs["RR"])

        if qrs["qt_duration"][chan] is not None:
            qt_intervals.append(qrs["qt_duration"][chan])

    print("Комплексы: {}".format(len(metadata)))
    print("Зубцы: {}".format(wcount))

    print("Типы комплексов:")
    print(classes)
    print("Автоклассы:")
    print(classids)

    print("ST-интервалы: {}, в среднем {} мс".format(
        len(st_intervals), np.mean(st_intervals)
    ))

    print("QT-интервалы: {}, в среднем {} мс".format(
        len(qt_intervals), 1000 * np.mean(qt_intervals)
    ))

    print("RR-интервалы: {}, в среднем {} мс".format(
        len(rr_intervals), 1000 * np.mean(rr_intervals)
    ))


def show_qrs(filename, chan, lim):
    sig, header = ecgread(filename)
    fs = header["fs"]

    sig = fix_baseline(
        sig,
        fs=fs,
        bias_window_ms=1500
    )

    if fs != 250:
        print("Warning! fs={}".format(fs))

    if lim:
        lim = min(lim, sig.shape[0])
    else:
        lim = sig.shape[0]

    metadata, pant = qrs_detection(
        sig[:lim,:3],
        fs=header["fs"]
    )

    approx, ders = ddwt(sig[:lim, chan], num_scales=5)

    fig, axarr = plt.subplots(len(ders)+1, 1, sharex="col")

    axarr[len(ders)].plot(pant)

    for i, s in enumerate(ders):
        axarr[i].plot(s,"b")
        axarr[i].plot(approx[i],"g")

        for ncycle, qrs in enumerate(metadata):
            lb = int(qrs["qrs_start"] * fs)
            rb = int(qrs["qrs_end"] * fs)
            axarr[i].plot(np.arange(lb, rb), s[lb:rb], "r")

        axarr[i].grid()

    #ms_step = 1000
    #locs = np.arange(0, len(ders[0]), ms_to_samples(ms_step, fs))
    #labs = [time.strftime("%M:%S", time.gmtime(x * ms_step)) for x in locs]
    #axarr[-1].set_xticks(locs, labs)
    #ms_step = 100
    #locs = np.arange(0, len(ders[0]), ms_to_samples(ms_step, fs))
    #axarr[-1].set_xticks(locs, minor=True)

    print("Look at the plots")
    plt.show()


def show_waves(filename, chan, smp_from=0, smp_to=0, draw=False):
    sig, header = ecgread(filename)

    print(sig.shape)

    fs = header["fs"]
    if fs != 250:
        print("Warning! fs={}".format(fs))

    sig = fix_baseline(
        sig,
        fs=fs,
        bias_window_ms=1500
    )

    if smp_to:
        smp_to = min(smp_to, sig.shape[0])
    else:
        smp_to = sig.shape[0]
    s = sig[smp_from:smp_to, chan]

    metadata, foo = qrs_detection(
        sig[smp_from:smp_to,:],
        fs=header["fs"]
    )

    find_points(
        sig[smp_from:smp_to, :],
        fs=header["fs"],
        bias=header["baseline"],
        metadata=metadata
    )

    metadata_postprocessing(
        metadata,
        sig[smp_from:smp_to, :],
        header
    )

    if draw:
        plt.plot(s, "b")

        # Цвета для раскрашивания зубцов на графике
        pt_keys = {
            "q_pos": "g",
            "r_pos": "r",
            "r2_pos": "m",
            "s_pos": "b",
            "s2_pos": "c",
            "p_pos": "y",
            "t_pos": "k",
            "t_start": "k",
            "t_end": "k"}

        for ncycle, qrs in enumerate(metadata):

            lb = int(qrs["qrs_start"] * fs)
            rb = int(qrs["qrs_end"] * fs)
            iso = qrs["isolevel"][chan] + header["baseline"][chan]
            plt.plot([lb, rb], [iso]*2, "g:")

            for k in pt_keys:
                point = qrs[k][chan]
                if point is not None:
                    plt.scatter(point, s[point], c=pt_keys[k])

            if qrs["qrs_center"] > 120:
                break

            lb = qrs["t_start"][chan]
            rb = qrs["t_end"][chan]
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

    print("Общая длительность {} c".format(len(s)/fs))

    chss = [x["heartrate"] for x in metadata if x["heartrate"] is not None]
    if chss:
        print("ЧСС: мин. {:.2f}, макс. {:.2f}".format(min(chss), max(chss)))

    missing_hrt = [i for i,x in enumerate(metadata) if x["heartrate"] is None]
    if missing_hrt:
        print("Heartrate missing in beats\n{}".format(missing_hrt))

    #plt.xlim((200,700))

    ry = define_rythm(metadata)
    rtm = {}
    for r in ry:
        desc = r["desc"]
        rtm[desc] = rtm.get(desc,0) + 1

    print("Ритмы:")
    print(rtm)

    print("Ишемия...")
    m = define_ishemia_episodes(sig[smp_from:smp_to, :], header, metadata)
    if m:
        print(m)
        print("Число эпизодов: {}".format(len(m)))

    #for x in r:
    #    print("ритм {}: {} с".format(x["desc"], x["end"] - x["start"]))

    print_summary(metadata, chan)

    if draw:
        #f, ax = plt.subplots(1,2)
        #show_qt_hist(ax[0], metadata, "qt_duration")
        #show_qt_hist(ax[1], metadata, "qtc_duration")

        plt.grid(True)
        plt.show()


def show_qt_hist(ax, metadata, key):
    numch = len(metadata[0]["r_pos"])
    for chan in range(numch):
        hdqt = calculate_histogram(metadata, key, channel=chan, bins=10,
                                   censoring=False)
        x = [hdqt[0]["bin_left"]*1000]
        y = [0]
        for bin in hdqt:
            x.append(bin["bin_right"]*1000)
            y.append(bin["percent"])

        x.append(hdqt[-1]["bin_right"] * 1000)
        y.append(0)

        ax.step(x,y)
        ax.set_title("{}".format(sum(x["count"] for x in hdqt)))


def main():
    # Rh2022 = qr, noise
    # Rh2021 - Rs, extracyc
    # Rh2024 - p(q)RsT
    # Rh2025 = rs
    # Rh2010 - дрейф, шум, артефакты

    #filename = "/Users/arseniy/SERDECH/data/PHYSIONET/I16"
    #filename = "testI59.ecg"
    filename = "TestFromDcm.ecg"
    #filename = "/Users/arseniy/SERDECH/data/ROXMINE/Rh2021"

    if not filename.endswith(".ecg") and not os.path.isfile(filename + ".hea"):
        print("Файл не найден")
        return

    #show_filter_responses()

    #show_raw(
    #    filename,
    #    smp_from=0,
    #    smp_to=2000
    #)

    #show_decomposition(
    #    filename,
    #    chan=1,
    #    smp_from=0,
    #    smp_to=20000
    #)

    #show_qrs(
    #    filename,
    #    chan=6,
    #    lim=20000
    #)

    show_waves(
        filename,
        chan=1,#common_signal_names.index("V1"),
        smp_from=0,
        smp_to=0,
        draw=False
    )


if __name__ == "__main__":
    main()