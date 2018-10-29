#!/usr/bin/env python
#coding: utf-8


from scipy.signal import convolve, hann
import numpy as np
from scipy.fftpack import fft, ifft
from util import signal_channels

"""
    Расчет усредненного амплитудного спектра
"""


def mean_spectrum(x, aperture=1024, log_output=True):

    N = len(x)
    acc = None
    ham = np.array(hann(aperture))

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


def build_comb_filter(fs, n, att, base_freq=50.0, q=5.0):
    """
    Построение АЧХ гребенчатого режекторного фильтра
    :param fs: частота дискретизации
    :param n: число точек в спектре
    :param att: ослабление гармоник 0 - 1
    :param base_freq: частота первой гармоники
    :param q: ширина полосы задержания
    :return: f_grid, response
    """

    att = min(1.0, max(att, 0.0))
    response = np.ones(n)
    f_grid = np.arange(0.0, fs, float(fs)/n)

    for i, f in enumerate(f_grid):
        real_f = f if f <= fs/2 else fs-f
        for harm in np.arange(base_freq, fs/2, base_freq):
            d = (1.0 - att) * np.exp(-((real_f-harm)/q)**2)
            response[i] = min(response[i], 1.0 - d)

    return f_grid, response


def mains_filter(sig, fs, bias, mains, attenuation, aperture):
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

    result = np.zeros(sig.shape)

    for chan, x in signal_channels(sig):

        y = np.zeros(len(x))
        hamwnd = np.array(hann(aperture))
        step = int(aperture / 2)

        ham_left = hamwnd.copy()
        ham_left[:step+1] = np.max(hamwnd)

        n2 = len(x) - aperture
        for n1 in range(0, n2, step):
            # для ослабления краевых эффектов берем несимметричное окно в
            # начале
            wnd = ham_left if n1 == 0 else hamwnd
            # комплексный спектр с учетом окна и смещения
            xf = fft(wnd * (np.array(x[n1:n1 + aperture]) - bias[chan]))
            # отфильтрованный сигнал в окне
            yt = np.real(ifft(xf * fft_response))

            y[n1:n1 + aperture] += yt

        result[:,chan] = y + bias[chan]

    return result


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
    h = hann(int(samples_per_ms * bias_window_ms))
    h = h / sum(h)

    result = np.zeros(sig.shape)

    for chan, x in signal_channels(sig):
        bias = np.mean(x)
        # огибающая (фон) вычисляется путем свертки со сглаживающей
        # апертурой и затем вычитается из входного сигнала
        bks = convolve(x - bias, h, mode="same")
        result[:,chan] = x - bks

    return result
