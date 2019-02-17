# coding: utf-8

from scipy.fftpack import fft
import numpy as np
from config import config
from metadata import is_artifact, is_pvc

from matplotlib import pyplot as plt


def resample_rhythm(metadata, dx):
    """resample_rhythm

    :param metadata:
    :param dx:
    :return: Пододвинутые значения RR в мс
    """
    sinc_wnd = 5

    buf = []
    rdpos = 0

    resamp = []

    for pos in np.arange(0, metadata[-1]["qrs_center"]+dx, dx):

        m = buf[0][0] if buf else -sinc_wnd
        while m < pos-sinc_wnd and buf:
            buf.pop(0)
            m = buf[0][0]

        m = metadata[rdpos]["qrs_center"]
        while m < pos+sinc_wnd:

            if not (is_artifact(metadata[rdpos]) or is_pvc(metadata[rdpos])):
                buf.append((m, 1000.0 * metadata[rdpos]["RR"]))

            rdpos += 1
            if rdpos >= len(metadata):
                break

            m = metadata[rdpos]["qrs_center"]

        if rdpos >= len(metadata):
            break

        s = 0
        for x, y in buf:
            s += np.sinc(x-pos)
        resamp.append(s)

    plt.plot()

    return np.array(resamp)


def plot_rhythm_resamp(metadata, dx):

    x0 = []
    y0 = []
    for x in metadata:
        if not (is_artifact(x) or is_pvc(x)):
            x0.append(x["qrs_center"])
            y0.append(1000.0 * x["RR"])

    y1 = resample_rhythm(metadata, dx=dx)
    x1 = np.arange(0, dx*len(y1), dx)

    plt.stem(x0, y0)
    plt.plot(x1, y1, "r")
    plt.show()


def rhythm_spectrum(metadata, **kwargs):
    """Расчет частотных составляющих сердечного ритма

    :param metadata: список метаданных (только чтение). Требуемые ключи:
        - qrs_center
        - RR
    :param kwargs:
        - freq_bands - список границ частотных диапазонов (Гц)
        - sampling - желаемая частота дискретизации сигнала ритма.
    :return:
        - results - словарь с оцениваемыми параметрами
        - retn - относительная мощность спектра по заданных диапазонам
        - fp - частоты для построения амплитудного спектра
        - sp - значения для построения амплитудного спектра []
    """

    bands = kwargs.get(
        "freq_bands",
        config.RSVAR["freq_bands"]
    )

    fs = kwargs.get(
        "sampling",
        config.RSVAR["sampling"]
    )

    r = resample_rhythm(metadata, dx=1/fs)

    n = len(r)
    sp = np.sqrt(1.0 / n) * np.abs(fft(r))

    fp = np.arange(0.0, fs, float(fs) / n)

    df = fp[1]
    ret = []
    totalpw = 0.0

    for i in range(1, len(bands)):

        f_lo = int(np.ceil(bands[i-1]/df))
        f_hi = int(np.ceil(bands[i]/df))
        pw = sum(np.array(sp[f_lo:f_hi])**2)
        totalpw += pw
        ret.append((bands[i-1], bands[i], pw))

    # нормировка по общей мощности
    retn = []
    for f1, f2, pw in ret:
        retn.append((f1, f2, pw/totalpw))

    tp = sum(sp*sp)  # общая мощность
    ulf = ret[0][2]  # мощность в диапазоне < 0,003 Гц
    vlf = ret[1][2]  # мощность в диапазоне 0,003 - 0,04 Гц
    tpm = tp - ulf - vlf
    lf = ret[2][2]  # мощность в диапазоне 0,04 - 0,15 Гц
    hf = ret[3][2]  # мощность в диапазоне 0,15 - 0,4 Гц

    results = {
        "TP": tp,
        "ULF": ulf,
        "VLF": vlf,
        "LF": lf,
        "HF": hf,
        "LFn": 100.0 * lf / tpm,
        "HFn": 100.0 * hf / tpm,
        "LFHF": lf / hf
    }

    m = int(n/2)
    return results, retn, fp[:m], sp[:m]
