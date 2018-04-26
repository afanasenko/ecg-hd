#!/usr/bin/env python
#coding: utf-8


from scipy import signal
import numpy as np
from scipy.fftpack import fft, ifft


def signal_channels(signal):
    if len(signal.shape) == 1:
        yield signal
    else:
        for channel in range(signal.shape[0]):
            yield signal[channel, :]


"""
    Поиск границ кардиоциклов
"""


def splitter(x):

    # убираем короткие выбросы
    h = signal.hann(100)
    xx = signal.convolve(x, h/sum(h), mode="same")

    mx = signal.argrelmax(xx, 0, 250)

    print("average period is {} samples".format(np.mean(np.diff(mx[0]))))

    return mx[0]

"""
    Вычисление спектра, усредненного по кардиоциклам
"""


def synchro_spectrum(xlo, xhi, log_output=False):

    # пики определяются по НЧ каналу
    mx = splitter(xlo)

    per = np.mean(np.diff(mx))

    # делаем небольшое перекрытие соседних окон
    hw = int(0.75 * (per-1))

    wnd = signal.hann(2*hw+1)

    agg_hi = np.zeros(hw+1)
    agg_lo = np.zeros(hw+1)
    aggc = 0


    for mm in mx:
        if mm >= hw and mm+hw < len(xhi):
            aggc += 1

            s = xhi[mm-hw : mm+hw+1]
            sp = fft(s * wnd)
            sp = np.abs(sp[0:hw+1])
            agg_hi += sp

            s = xlo[mm - hw: mm + hw + 1]
            sp = fft(s * wnd)
            sp = np.abs(sp[0:hw + 1])
            agg_lo += sp

    agg_lo /= aggc
    agg_hi /= aggc

    if log_output:
        elo = 20000#np.sqrt(sum(np.power(agg_lo, 2)))
        agg_lo = 20.0 * np.log10(agg_lo / elo)

        ehi = 10000#np.sqrt(sum(np.power(agg_hi, 2)))
        agg_hi = 20.0 * np.log10(agg_hi / ehi)

    return agg_lo, agg_hi


"""
    Вычисление усредненного кардиоцикла
"""


def synchro_cycle(xlo, xhi):

    # пики определяются по НЧ каналу
    mx = splitter(xlo)

    per = np.mean(np.diff(mx))

    # делаем небольшое перекрытие соседних окон
    hw = int(0.75 * (per-1))

    wnd = signal.hann(2*hw+1)

    agg_hi = np.zeros(2*hw+1)
    agg_lo = np.zeros(2*hw+1)
    aggc = 0


    for mm in mx:
        if mm >= hw and mm+hw < len(xhi):
            aggc += 1

            s = xhi[mm-hw : mm+hw+1]
            agg_hi += s

            s = xlo[mm - hw: mm + hw + 1]
            agg_lo += s

    print("signal was aggregated across {} cycles".format(aggc))

    return agg_lo / aggc, agg_hi / aggc

"""
    Вычисление спектра, усредненного по кардиоциклам
"""


def synchro_spectrum_r(sig, peaks, halfw, log_output=False):

    halfw = int(halfw)
    wnd = signal.hann(2*halfw+1)

    agg = np.zeros(halfw+1)
    aggc = 0


    for mm in peaks:
        if mm >= halfw and mm+halfw < len(sig):
            aggc += 1

            s = sig[mm-halfw : mm+halfw+1]
            sp = fft(s * wnd)
            sp = np.abs(sp[0:halfw+1])
            agg += sp

    agg /= aggc

    if log_output:
        elo = 20000#np.sqrt(sum(np.power(agg_lo, 2)))
        agg = 20.0 * np.log10(agg / elo)

    return agg


def plain_spectrum_r(sig, peaks, halfw, log_output=False):

    wnd = signal.hann(2*halfw+1)

    agg = np.zeros(halfw+1)
    aggc = 0


    for mm in peaks:
        if mm >= halfw and mm+halfw < len(sig):
            aggc += 1

            s = sig[mm-halfw : mm+halfw+1]
            sp = fft(s * wnd)
            sp = np.abs(sp[0:halfw+1])
            agg += sp

    agg /= aggc

    if log_output:
        elo = 20000#np.sqrt(sum(np.power(agg_lo, 2)))
        agg = 20.0 * np.log10(agg / elo)

    return agg

"""
    Расчет усредненного амплитудного спектра
"""


def mean_spectrum(x, aperture=1024, log_output=True):

    N = len(x)
    acc = None
    ham = np.array(signal.hann(aperture))

    n1 = 0
    happ = int(aperture/2)
    K = 0
    while n1 < N - aperture:
        xs = np.array(x[n1:n1 + aperture])
        yf = fft(xs * ham)
        yf = np.abs(yf[0:happ])
        K += 1

        if acc is not None:
            acc += yf
        else:
            acc = yf

        n1 += happ

    acc /= K

    if log_output:
        elo = 20000#np.sqrt(sum(np.power(agg_lo, 2)))
        acc = 20.0 * np.log10(acc / elo)

    return acc

"""
    Предварительная обработка сигнала, включающая
    вычитание фона
    выделение R-зубцов
"""


def preprocess(signals):
    signals["unbias"] = []
    signals["hif"] = []
    signals["zub"] = []

    for i, sig in enumerate(signals["rawdata"]):
        # убираем фон
        h = signal.hann(1000)
        xx = sig - signal.convolve(sig, h / sum(h), mode="same")
        signals["unbias"].append(xx)

        h = signal.hann(80)
        yy = xx - signal.convolve(xx, h / sum(h), mode="same")
        #signals["hif"].append(xx.copy()/5) # 20151012
        signals["hif"].append(yy)

        mx = signal.argrelmax(yy, 0, 250)

        pks = []
        sigsize = len(signals["unbias"][-1])
        for x in mx[0]:
            if signals["hif"][-1][x] > 30:
                # уточняем по первичному сигналу
                n1 = max(0, x - 45)
                n2 = min(sigsize, x + 45)
                delta = np.argmax(signals["unbias"][-1][n1:n2])
                pks.append(n1 + delta)

        signals["zub"].append(pks)


def build_comb_filter(fs, n, att, base_freq=50.0, q=5.0):

    att = min(1.0, max(att, 0.0))
    response = np.ones(n)
    f_grid = np.arange(0.0, fs, float(fs)/n)

    for i, f in enumerate(f_grid):
        real_f = f if f <= fs/2 else fs-f
        for harm in np.arange(base_freq, fs/2, base_freq):
            d = (1.0 - att) * np.exp(-((real_f-harm)/q)**2)
            response[i] = min(response[i], 1.0 - d)

    return f_grid, response


def mains_filter(sig, fs, mains, attenuation, aperture):
    """
    Подавление гармоник частоты электрической сети
    :param sig:
    :param fs: частота дискретизации в Гц
    :param mains: частота сети
    :param attenuation: коэффициент ослабления гармоник (0 - полное подавление)
    :param aperture: апертура БПФ
    :return:
    """

    f_grid, fft_response = build_comb_filter(
        fs=fs,
        n=aperture,
        att=attenuation,
        base_freq=mains,
        q=mains*0.03
    )

    result = []

    for x in signal_channels(sig):

        y = np.zeros(len(x))
        hamwnd = np.array(signal.hann(aperture))
        step = int(aperture / 2)

        n2 = len(x) - aperture
        for n1 in range(0, n2, step):
            # комплексный спектр с учетом окна
            xf = fft(hamwnd * np.array(x[n1:n1 + aperture]))
            # отфильтрованный сигнал в окне
            yt = np.real(ifft(xf * fft_response))

            y[n1:n1 + aperture] += yt

        result.append(y)

    return np.array(result)


def fix_baseline(sig, fs, bias_window_ms):
    """
    fix_baseline выравнивает базовую линию
    :param sig: numpy array - отсчеты сигнала
    :param fs: частота дискретизации, Гц
    :param bias_window_ms: ширина окна для подаввления фона (мс)
    :return: сигнал с подавленным фоном
    """

    samples_per_ms = float(fs)/1000

    # косинусоидальная сглаживающая апертура шириной bias_window_ms
    h = signal.hann(int(samples_per_ms * bias_window_ms))
    h = h / sum(h)

    result = []

    for x in signal_channels(sig):
        bias = np.mean(x)
        # огибающая (фон) вычисляется путем свертки со сглаживающей апертурой и затем вычитается из входного сигнала
        bks = signal.convolve(x - bias, h, mode="same")
        result.append(x - bks)

    return np.array(result)
