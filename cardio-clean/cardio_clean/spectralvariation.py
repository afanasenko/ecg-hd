# coding: utf-8

from scipy.fftpack import fft
import numpy as np
from config import config
from metadata import is_artifact, is_pvc

from matplotlib import pyplot as plt


def resample_rhythm(metadata, dx):

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
                buf.append((m, 1.0 / metadata[rdpos]["RR"]))

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
            y0.append(1.0 / x["RR"])

    y1 = resample_rhythm(metadata, dx=dx)
    x1 = np.arange(0, dx*len(y1), dx)

    plt.stem(x0, y0)
    plt.plot(x1, y1, "r")
    plt.show()


def rhythm_spectrum(metadata, bands=(0.0, 0.003, 0.04, 0.15, 0.4), fs=1.0):

    plot_rhythm_resamp(metadata, dx=1/fs)

    r = resample_rhythm(metadata, dx=1/fs)

    n = len(r)
    m = int(n/2)
    sp = np.abs(fft(r))
    sp[0] = 0
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

    return retn, fp[:m], sp[:m]
