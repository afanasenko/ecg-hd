# coding: utf-8

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import lfilter, hann
from scipy.fftpack import fft

from graphresults import ecgread
from wavdetect import ddwt


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


def zcfind(x):
    return []


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


def zcfind(x):

    zc = []

    for i in range(1,len(x)):
        if x[i-1]*x[i] < 0:
            w1 = float(abs(x[i-1]))
            w2 = float(abs(x[i]))
            zc.append(i-0.5)

    return np.array(zc)


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
        for x in ders:
            print(zcfind(x) - mid)


def show_decomposition():

    sig, hdr = ecgread("/Users/arseniy/SERDECH/data/ROXMINE/Rh2022")
    ders = ddwt(sig[:1200,0], num_scales=5)
    multiplot(ders)

if __name__ == "__main__":

    #show_filter_responses()
    show_decomposition()

