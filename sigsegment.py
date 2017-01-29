#!/usr/bin/env python
#coding: utf-8

from scipy import signal
import numpy as np
import wfdb
import matplotlib.pyplot as plt
import sys

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
        lfsignal = x - signal.convolve(x, h, mode="same")
    else:
        lfsignal = np.array(x, 'float')

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
    
    search_window = samples_per_ms * 10  # 10 миллисекунд

    for pos in extrema[0]:
        if hfsignal[pos] > peak_threshold:
            # уточняем максимум по первичному сигналу, просматривая окрестности текущего отсчета
            n1 = int(max(0, pos - search_window))
            n2 = int(min(sigsize, pos + search_window))
            delta = np.argmax(lfsignal[n1:n2])
            pks.append(n1 + delta)

    # результат можно преобразовать в миллисекунды по формуле 1000 * pks / fs
    return np.array(pks) 

def main():
    # загрузка данных
    recordname = sys.argv[1] # The name of the WFDB record to be read (without any file extensions)
    sampto = int(sys.argv[2]) # The final sample number to read for each channel
    peak_interval_ms = int(sys.argv[3])
    sig, fields=wfdb.rdsamp(recordname, sampto=sampto)
    x = sig[:,0]
    fs = fields["fs"] # sampling frequency (in samples per second per signal)

    pks = extract_short_peaks(x, fs, peak_interval_ms=peak_interval_ms)
    A = [x[i] for i in pks]
    I = [(pks[i+1] - pks[i])/1000 for i in range(len(pks)-1)]

    plt.style.use("ggplot")
    t=np.array(range(0,sig.shape[0]))/fs
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 12
    fig_size[1] = 6
    plt.rcParams["figure.figsize"] = fig_size
    plt.plot(t, x)
    plt.plot(pks/fs, A, "b+", markersize=6)
    #plt.plot(pks/fs, np.zeros(len(pks)), "gd", markersize=6)
    plt.title("MIT-BIH Arrhythmia record")
    plt.xlabel("time/s")
    plt.ylabel(fields["units"][0])
    plt.show()
if __name__ == "__main__":
    main()