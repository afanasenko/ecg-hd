#!/usr/bin/env python
#coding: utf-8

from scipy import signal
import numpy as np
import wfdb
import matplotlib.pyplot as plt
import sys


def erode(x,n):
    """
    Эрозия обычная
    :param x:
    :param n:
    :return:
    """

    r = max(1, np.floor(n/2))
    output = np.zeros(x.shape)

    for i in range(len(x)):
        i1 = max(0, i-n)
        i2 = min(len(x)-1, i+n)
        output[i] = min(x[i1:i2])

    return output


def conditional_dilate(x, mask, n=3):
    """
     Условное наращивание. Выходной сигнал не превосходит значений сигнала mask
     :param x:
     :param mask:
     :param n: ширина структурного элемента
     :return: (output, changed). Changed - флаг того, что сигнал изменился
    """

    r = int(max(1, np.floor(n/2)))
    output = np.zeros(x.shape)
    changed = False
    nc = 0

    for i in range(len(x)):
        i1 = max(0, i-r)
        i2 = min(len(x)-1, i+r)
        dil = max(x[i1:i2])
        # сравнение с маской
        dil = min(dil, mask[i])
        if dil != x[i]:
            changed = True
            nc += 1

        output[i] = dil

    #sys.stdout.write("{}, ".format(nc))
    return output, changed

def open_by_reconstruction(x, strel_size):
    """
     Морфологическое открытие через реконструкцию. "Срезает" положительные выбросы
     :param x:
     :param fs:
     :return:
    """

    marker = erode(x, strel_size)
    progress = True

    while progress:
        marker, progress = conditional_dilate(marker, mask=x, n=strel_size)

    return marker


def extract_peaks_morpho(x, fs, peak_length_ms=20):
    samples_per_ms = fs/1000
    elsize = max(3, round(samples_per_ms * peak_length_ms))
    return open_by_reconstruction(x, elsize)


def extract_short_peaks(x, fs, bias_window_ms=250, peak_length_ms=20, peak_interval_ms=500):
    """
    extract_short_peaks локализует всплески сигнала с продолжительностью меньше заданной
    :param x: numpy array - отсчеты сигнала
    :param fs: частота дискретизации, Гц
    :param bias_window_ms: ширина окна для подаввления фона (мс). 0 - фон не подавляется
    :param peak_length_ms: максимальная длительность всплеска
    :param peak_interval_ms: минимальный интервал времени между всплесками (для гашения помех)
    :return: np_array с номерами отсчетов сигнала, в которых найдены всплески
    """

    samples_per_ms = fs/1000
    #print(samples_per_ms * peak_interval_ms)

    if bias_window_ms:
        # косинусоидальная сглаживающая апертура шириной bias_window_ms
        h = signal.hann(samples_per_ms * bias_window_ms)
        h = h / sum(h)

        # огибающая (фон) вычисляется путем свертки со сглаживающей апертурой и затем вычитается из входного сигнала
        bks = signal.convolve(x, h, mode="same")
        lfsignal = x - bks
    else:
        lfsignal = np.array(x, 'float')
        bks = None

    # Готовим высокочастотный препарат, подчеркивающий короткие выбросы

    h = signal.hann(samples_per_ms * peak_length_ms)
    h = h / sum(h)

    hfsignal = lfsignal - signal.convolve(lfsignal, h, mode="same")
    #print(hfsignal)
    # по данному ВЧ препарату находим локальные максимумы, отстоящие друг от друга не менее, чем на peak_interval_ms
    extrema = signal.argrelmax(hfsignal, 0, int(samples_per_ms * peak_interval_ms))

    pks = []
    sigsize = len(x)

    # А дальше начинается самая колхозная часть: мы нашли максимумы в "искаженном" сигнале (hfsignal),
    # а они могут быть сдвинуты относительно пиков исходного сигнала x.
    # Поэтому мы будем просматривать окрестности каждого пика.
    # Кроме того, сигнал hfsignal по определению более шумный, и локальные максимумы в нем могут быть ложными.
    # Поэтому мы введем порог на значение сигнала в пике.

    peak_threshold = 0   # размерность не определена и вообще это самый подлый параметр, от него надо избавляться.
    # для хороших сигналов можно попробовать порог 0, он универсальный
    
    search_window = samples_per_ms * 100  # 10 миллисекунд

    for pos in extrema[0]:
        if hfsignal[pos] > peak_threshold:
            # уточняем максимум по первичному сигналу, просматривая окрестности текущего отсчета
            n1 = int(max(0, pos - search_window))
            n2 = int(min(sigsize, pos + search_window))
            delta = np.argmax(lfsignal[n1:n2])
            pks.append(n1 + delta)

    # результат можно преобразовать в миллисекунды по формуле 1000 * pks / fs
    return np.array(pks), bks, hfsignal


def main_new(recordname, chan, show):
    print("processing " + recordname)
    # загрузка данных
    outname = recordname.split("/")[-1] + "_features.tsv"

    sampto = 3000 # The final sample number to read for each channel
    peak_length = 20
    peak_interval_ms = 690
    bias_ms = 1.9*peak_interval_ms
    #sig, fields=wfdb.rdsamp(recordname, sampto=sampto)
    sig, fields=wfdb.rdsamp(recordname, sampto=450000)
    x = sig[:, chan]
    x -= np.mean(x)
    fs = fields["fs"] # sampling frequency (in samples per second per signal)
    print("Sampling frequency: {} Hz".format(fs))
    print("Duration: {} s".format(len(x) / fs))

    basel = extract_peaks_morpho(x, fs, peak_length_ms=3*peak_length)

    # немного сглаживаем биномиальным фильтром
    binfilt = [0.25, 0.5, 0.25]

    prefilt = signal.convolve(x - basel, binfilt, mode="same")

    rpeaks, bks, hfs = extract_short_peaks(x, fs, bias_window_ms=bias_ms, peak_interval_ms=peak_interval_ms)


    dist, bins = np.histogram(prefilt, 25)
    #axarr[1].plot(bins[:-1], dist, "r")
    q = np.percentile(prefilt, [87, 90, 95])
    print("\n")
    print(q)

    cycvalid = 0

    cdata = []

    tpeaks = []
    tamp = []
    ramp = []

    for ncyc in range(len(rpeaks)-1):
        cyc_begin = rpeaks[ncyc]
        cyc_end = min(len(prefilt), rpeaks[ncyc+1])
        state = 0

        tw = [0, 0]
        rw = [cyc_begin, 0]
        pw = [0, 0]
        tt = 0.095#q[0]
        ttm = 0.2*tt
        ta = 0
        ra = 0

        for i in range(cyc_begin, cyc_end):
            # R-wave
            if state == 0:
                if prefilt[i] < tt:
                    state = 1
                    rw[1] = i
                    ra = prefilt[cyc_begin]
            # R-T segment
            elif state == 1:
                if prefilt[i] > tt:
                    state = 2

                    j = i
                    while j>rw[1] and prefilt[j]>ttm:
                        j -= 1
                    tw[0] = j

            # T-wave
            elif state == 2:
                if prefilt[i] < tt:
                    state = 3
                    tw[1] = i
                    ta = max(prefilt[tw[0]:tw[1]])
                    tp = tw[0] + np.argmax(prefilt[tw[0]:tw[1]])
                    tpeaks.append(tp)

            # T-P segment
            elif state == 3:
                if prefilt[i] > tt:
                    state = 4
                    pw[0] = i
                # задний фронт T-зубца
                j = tw[1]
                while j < i and prefilt[j] > ttm:
                    j += 1
                tw[1] = j

            # P-wave
            elif state == 4:
                if prefilt[i] < tt:
                    state = 5
                    pw[1] = i
            # P-Q segment
            elif state == 5:
                if prefilt[i] > tt:
                    state = 6

        if state == 6:
            cycvalid += 1
            #print("cycle {:04}: state {}, tdur {}".format(ncyc, state, tw[1] - tw[0]))
            cdata.append([
                ra,
                (rw[1] - rw[0]) / fs,
                ta,
                (tw[1]-tw[0])/fs
            ])

            ramp.append(ra)
            tamp.append(ta)

    with open(outname, "w") as fp:
        for samp in cdata:
            fp.write("\t".join([str(x) for x in samp]) + "\n")

    print("{} of {} T-waves recognized".format(cycvalid, len(rpeaks)))

    if show:
        #plt.style.use("ggplot")
        t = np.arange(0, sampto) / fs
        fig_size = plt.rcParams["figure.figsize"]
        plt.rcParams["figure.facecolor"] = "white"
        fig_size[0] = 12
        fig_size[1] = 6
        plt.rcParams["figure.figsize"] = fig_size
        fig, axarr = plt.subplots(1, 1)
        axarr.plot(t, x[:sampto], 'r')
        plt.hold(True)

        rp_show = np.array([x for x in rpeaks if x < sampto])
        rp_val = [x[i] for i in rp_show]

        tp_show = np.array([x for x in tpeaks if x < sampto])
        tp_val = [x[i] for i in tp_show]

        axarr.plot(rp_show / fs, rp_val, "k*", markersize=6)
        axarr.plot(tp_show / fs, tp_val, "k+", markersize=6, linewidth=2)
        #axarr[0].plot(t, prefilt[:sampto], "g", alpha=0.5)
        #axarr[0].set_title("ECG segmentation: {}".format(recordname.split('/')[-1]))
        axarr.set_xlabel("время, с")
        axarr.set_ylabel("напряжение, мВ")

        #axarr[1].scatter(tamp, ramp)

        print("look at the plots")
        plt.show()


def main_showstages(recordname, chan):
    print("processing " + recordname)


    sampto = 3000 # The final sample number to read for each channel
    peak_length = 20
    peak_interval_ms = 690
    bias_ms = 1.9*peak_interval_ms
    peak_length_ms = 15
    #sig, fields=wfdb.rdsamp(recordname, sampto=sampto)
    sig, fields=wfdb.rdsamp(recordname, sampto=450000)
    x = sig[:, chan]
    x -= np.mean(x)
    fs = fields["fs"] # sampling frequency (in samples per second per signal)

    #basel = extract_peaks_morpho(x, fs, peak_length_ms=3*peak_length)

    # немного сглаживаем биномиальным фильтром
    #binfilt = [0.25, 0.5, 0.25]

    #prefilt = signal.convolve(x - basel, binfilt, mode="same")

    rpeaks, bks, hfs = extract_short_peaks(
        x,
        fs,
        bias_window_ms=bias_ms,
        peak_length_ms=peak_length_ms,
        peak_interval_ms=peak_interval_ms
    )


    dist, bins = np.histogram(x, 25)
    #axarr[1].plot(bins[:-1], dist, "r")
    q = np.percentile(x, [87, 90, 95])
    print("\n")
    print(q)

    #plt.style.use("ggplot")
    t = np.arange(0, sampto) / fs
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 12
    fig_size[1] = 6
    plt.rcParams["figure.figsize"] = fig_size
    plt.rcParams["figure.facecolor"] = "white"
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212, sharex=ax1)

    ax1.plot(t, x[sampto:2*sampto], 'r')
    plt.hold(True)

    #rp_show = np.array([x for x in rpeaks if x < sampto])
    #rp_val = [x[i] for i in rp_show]

    #tp_show = np.array([x for x in tpeaks if x < sampto])
    #tp_val = [x[i] for i in tp_show]

    #axarr.plot(rp_show / fs, rp_val, "b+", markersize=6)
    #axarr.plot(tp_show / fs, tp_val, "m+", markersize=6)
    ax1.plot(t, bks[sampto:2*sampto], "k", linewidth=2)
    ax2.plot(t, hfs[sampto:2*sampto], "k")
    #axarr.set_title("ECG segmentation: {}".format(recordname.split('/')[-1]))
    ax1.set_xlabel("время, с")
    ax1.set_ylabel("напряжение, мВ")
    ax1.legend(["входной сигнал", "низкочастотный фон"])
    ax2.set_xlabel("время, с")
    ax2.set_ylabel("напряжение, мВ")

    print("look at the plots")
    plt.show()

if __name__ == "__main__":
    main_new('./ecg_data/e0118', 0, True)

    #main_showstages('./ecg_data/e0602', 0)