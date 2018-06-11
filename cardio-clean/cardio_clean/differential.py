# coding: utf-8

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import lfilter, hann
from scipy.fftpack import fft

from graphresults import ecgread
from wavdetect import ddwt, find_points, zcfind
from qrsdetect import qrs_detection

def dummy_shift(x, n):
    return np.concatenate((x[n:], np.ones(n)*x[-1]))


def ddgau(sigma, fs, ksigma=3):

    dt = 1.0 / fs
    halfw = sigma*ksigma

    x = np.arange(-halfw, halfw, dt)
    y = np.zeros(len(x))
    sumsq = 0.0

    for i, t in enumerate(x):
        tt = t*t/(2.0*sigma*sigma)
        y[i] = (1-(t/sigma)**2)*np.exp(-tt)
        sumsq += y[i]*y[i]

    return x, y/sumsq


def dgau(sigma, fs, ksigma=3):

    dt = 1.0 / fs
    halfw = sigma*ksigma

    x = np.arange(-halfw, halfw, dt)
    y = np.zeros(len(x))
    sumsq = 0.0

    for i, t in enumerate(x):
        tt = t/sigma
        y[i] = -t*np.exp(-tt*tt)
        sumsq += y[i]*y[i]

    return x, y/sumsq


def gaussogram(x):

    sigmas = np.array([
        0.009,
        0.015,
        0.03,
        0.06
    ])

    fig, axarr = plt.subplots(len(sigmas)+1, 1, sharex=True)

    filterbank = []
    for i, s in enumerate(sigmas):
        t, h = dgau(s,250)
        filterbank.append(h)

    axarr[0].plot(x)

    for i, h in enumerate(filterbank):

        sm = lfilter(h, [1.0], x)
        sm = dummy_shift(sm, int(len(h)/2))

        axarr[i+1].plot(sm)
        axarr[i+1].grid()
    plt.show()


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
    ders = ddwt(x, num_scales=5)

    if spectral:
        multilinespec(ders[1:])
    else:
        multiplot(ders)
        for i, x in enumerate(ders):
            dt_theor = (2.0**i - 1)/2
            dt_exrep = zcfind(x) - mid
            print(dt_theor, dt_exrep)


def show_decomposition():

    sig, hdr = ecgread("/Users/arseniy/SERDECH/data/ROXMINE/Rh1001")
    ders = ddwt(sig[:1200,0], num_scales=3)
    multiplot(ders)


def show_waves():
    # Rh2022 = qr
    # Rh2021 - Rs, extracyc
    # Rh2025 = rs
    sig, header = ecgread("/Users/arseniy/SERDECH/data/ROXMINE/Rh2022")
    fs = header["fs"]
    if fs != 250:
        print("Warning! fs={}".format(fs))

    s = sig[:,0]

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

        if "rwav" in qrs:
            rwav = qrs["rwav"]
            plt.plot(np.arange(rwav[0],rwav[1]), s[rwav[0]:rwav[1]], "r")

        qp = qrs.get("qWavePosition", 0)
        if qp:
            plt.scatter(qp, s[qp])

        rp = qrs.get("rWavePosition", 0)
        if rp:
            plt.scatter(rp, s[rp])

        sp = rp = qrs.get("sWavePosition", 0)
        if sp:
            plt.scatter(sp, s[sp])

    plt.xlim((200,700))

    print(qrsTypes)

    plt.show()

if __name__ == "__main__":

    #show_filter_responses()
    show_decomposition()
    #show_waves()

