# coding: utf-8

""" Обнаружение сигналов каридиостимулятора
"""

import os
import numpy as np

from metadata import *
from scipy.signal import convolve, hann, argrelmax
from matplotlib import pyplot as plt

from cardio_clean.util import ecgread, signal_channels


def erode(x, n):
    """ Эрозия обычная
    :param x:
    :param n:
    :return:
    """

    r = int(max(1, np.floor(n / 2)))
    output = np.zeros(x.shape)

    for i in range(len(x)):
        i1 = max(0, i - r)
        i2 = min(len(x), i + r + 1)
        output[i] = min(x[i1:i2])

    return output


def conditional_dilate(x, mask, n=3):
    """ Условное наращивание. Выходной сигнал не превосходит значений
    сигнала mask
     :param x:
     :param mask:
     :param n: ширина структурного элемента
     :return: (output, changed). Changed - флаг того, что сигнал изменился
    """

    r = int(max(1, np.floor(n / 2)))
    output = np.zeros(x.shape)
    changed = False
    nc = 0

    for i in range(len(x)):
        i1 = max(0, i - r)
        i2 = min(len(x), i + r + 1)
        dil = max(x[i1:i2])
        # сравнение с маской
        dil = min(dil, mask[i])
        if dil != x[i]:
            changed = True
            nc += 1

        output[i] = dil

    return output, changed


def open_by_reconstruction(x, strel_size):
    """ Морфологическое открытие через реконструкцию. "Срезает"
    положительные выбросы
     :param x:
     :param fs:
     :return:
    """

    marker = erode(x, strel_size)
    progress = True

    while progress:
        marker, progress = conditional_dilate(marker, mask=x, n=strel_size)

    return marker


def hats(x, fs, peak_length_ms, strel_radius=None):
    """
     :param x:
     :param fs:
     :return:
    """

    if strel_radius is None:
        samples_per_ms = fs / 1000
        elsize = max(3, round(samples_per_ms * peak_length_ms))
    else:
        elsize = 1 + 2*strel_radius

    marker = erode(x, elsize)
    progress = True

    while progress:
        marker, progress = conditional_dilate(marker, mask=x, n=3)

    return x - marker


def extract_peaks_morpho(x, fs, peak_length_ms):
    samples_per_ms = fs / 1000
    elsize = max(3, round(samples_per_ms * peak_length_ms))

    pk = x - open_by_reconstruction(x, elsize)
    return pk


def extract_short_peaks(x, fs, bias_window_ms=250, peak_length_ms=20,
                        peak_interval_ms=500):
    """
    extract_short_peaks локализует всплески сигнала с продолжительностью меньше заданной
    :param x: numpy array - отсчеты сигнала
    :param fs: частота дискретизации, Гц
    :param bias_window_ms: ширина окна для подаввления фона (мс). 0 - фон не подавляется
    :param peak_length_ms: максимальная длительность всплеска
    :param peak_interval_ms: минимальный интервал времени между всплесками (для гашения помех)
    :return: np_array с номерами отсчетов сигнала, в которых найдены всплески
    """

    samples_per_ms = float(fs) / 1000
    # print(samples_per_ms * peak_interval_ms)

    if bias_window_ms:
        # косинусоидальная сглаживающая апертура шириной bias_window_ms
        h = hann(int(samples_per_ms * bias_window_ms))
        h = h / sum(h)

        # огибающая (фон) вычисляется путем свертки со сглаживающей апертурой и затем вычитается из входного сигнала
        bks = convolve(x, h, mode="same")
        lfsignal = x - bks
    else:
        lfsignal = np.array(x, 'float')
        bks = None

    # Готовим высокочастотный препарат, подчеркивающий короткие выбросы

    h = hann(int(samples_per_ms * peak_length_ms))
    h = h / sum(h)

    hfsignal = lfsignal - convolve(lfsignal, h, mode="same")
    # print(hfsignal)
    # по данному ВЧ препарату находим локальные максимумы, отстоящие друг от друга не менее, чем на peak_interval_ms
    extrema = argrelmax(hfsignal, 0,
                               int(samples_per_ms * peak_interval_ms))

    pks = []
    sigsize = len(x)

    # А дальше начинается самая колхозная часть: мы нашли максимумы в "искаженном" сигнале (hfsignal),
    # а они могут быть сдвинуты относительно пиков исходного сигнала x.
    # Поэтому мы будем просматривать окрестности каждого пика.
    # Кроме того, сигнал hfsignal по определению более шумный, и локальные максимумы в нем могут быть ложными.
    # Поэтому мы введем порог на значение сигнала в пике.

    peak_threshold = 0  # размерность не определена и вообще это самый подлый параметр, от него надо избавляться.
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
    return np.array(pks), bks, hfsignal, lfsignal


def demo_pm():
    #filename = "/Users/arseniy/SERDECH/data/PHYSIONET/102"
    filename = "/Users/arseniy/SERDECH/data/PHYSIONET/I59"

    chan = 0
    smp_from = 0
    smp_to = 20000

    if not filename.endswith(".ecg") and not os.path.isfile(filename + ".hea"):
        print("Файл не найден")
        return

    sig, header = ecgread(filename)
    fs = header["fs"]

    s = sig[smp_from:smp_to, chan]
    #s = np.abs(s)

    pk1 = hats(s, fs, peak_length_ms=2, strel_radius=1)
    pk2 = hats(s, fs, peak_length_ms=2, strel_radius=2)
    pk3 = hats(s, fs, peak_length_ms=2, strel_radius=3)

    t = np.arange(smp_from, smp_to) * 1.0 / fs

    plt.plot(t, s, t, pk1, t, pk2, t, pk3)

    plt.show()

if __name__ == "__main__":

    demo_pm()