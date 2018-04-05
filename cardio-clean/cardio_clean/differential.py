# coding: utf-8

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import lfilter

from graphresults import ecgread


def dummy_shift(x, n):
    return np.concatenate((np.ones(n)*x[-1], x[n:]))
    #return np.concatenate((x[n:], np.ones(n)*x[-1]))


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
        0.01,
        0.02,
        0.03,
        0.04
    ])

    fig, axarr = plt.subplots(len(sigmas)+1, 1, sharex=True)

    filterbank = []
    maxn = 0
    for i, s in enumerate(sigmas):
        t, h = dgau(s,250)
        filterbank.append(h)
        maxn = max(maxn, int(len(h)/2))
        print(maxn)

    axarr[0].plot(dummy_shift(x, maxn))

    for i, h in enumerate(filterbank):

        sm = lfilter(h, [1.0], x)
        sm = dummy_shift(sm, maxn - int(len(h)/2))

        axarr[i+1].plot(sm)
        axarr[i+1].grid()
    plt.show()


if __name__ == "__main__":
    #x, y = dgau(0.2, 250)
    #plt.plot(x, y)
    #plt.show()

    sig, hdr = ecgread("/Users/arseniy/SERDECH/data/ROXMINE/Rh2011")

    gaussogram(sig[:20000, 0])

