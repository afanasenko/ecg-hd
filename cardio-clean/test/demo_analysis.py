# coding: utf-8

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import lfilter, hann
from scipy.fftpack import fft

from cardio_clean.wavdetect import ddwt, find_points, zcfind
from cardio_clean.qrsdetect import qrs_detection
from demo_preprocessing import ecgread


def dummy_shift(x, n):
    return np.concatenate((x[n:], np.ones(n)*x[-1]))


def multiplot(siglist):
    fig, axarr = plt.subplots(len(siglist), 1, sharex=True)
    for i, s in enumerate(siglist):
        axarr[i].plot(s)
        axarr[i].grid()
    print("Look at the plots")
    plt.show()


def multispectrum(siglist):
    fig, axarr = plt.subplots(len(siglist), 1, sharex=True)
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


def show_waves():
    # Rh2022 = qr, noise
    # Rh2021 - Rs, extracyc
    # Rh2024 - p(q)RsT
    # Rh2025 = rs
    sig, header = ecgread("/Users/arseniy/SERDECH/data/ROXMINE/Rh1011")
    fs = header["fs"]
    if fs != 250:
        print("Warning! fs={}".format(fs))

    s = sig[:20000,0]

    metadata, debugdata = qrs_detection(
        sig[:,:],
        fs=header["fs"],
        bias=header["baseline"],
        gain=header["adc_gain"],
        minqrs_ms=20)

    newmeta = find_points(s, header["fs"], metadata)

    plt.plot(s, "b")

    qrsTypes = {}

    for qrs in newmeta:

        qrstype = qrs["qrsType"]
        qrsTypes[qrstype] = qrsTypes.get(qrstype, 0) + 1

        for k,v in qrs.get("waves", {}).items():
            c = v["center"]
            if c is not None:
                plt.scatter(c, s[c])

            if k in ("t", "p"):
                lb = v["start"]
                rb = v["end"]
                if lb is not None and rb is not None:
                    plt.plot(np.arange(lb, rb), s[lb:rb], "r")

    plt.xlim((200,700))

    print(qrsTypes)

    plt.show()

if __name__ == "__main__":

    #show_filter_responses()
    #show_decomposition()
    show_waves()

