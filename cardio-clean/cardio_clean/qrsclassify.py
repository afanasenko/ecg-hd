# coding: utf-8

import numpy as np
from scipy.signal import lfilter

from sigbind import signal_channels


def correlate(sig1, sig2):
    """
        Вычисление статического коэффициента корреляции многомерных сигналов
    :param sig1:
    :param sig2:
    :return:
    """

    samples = sig1.shape[0]

    sx = 0.0
    sy = 0.0
    sxy = 0.0
    for i in range(samples):
        sxy += np.dot(sig1[i,:], sig2[i,:])
        sx += np.dot(sig1[i,:], sig1[i,:])
        sy += np.dot(sig2[i,:], sig2[i,:])

    return sxy / (np.sqrt(sx) * np.sqrt(sy))


def newqrs():
    pass


def get_qrs_bounds(meta, fs):
    left = int(meta["qrs_start"] / fs)
    right = int(meta["qrs_end"] / fs)
    center = int(meta["r_wave_center"][0] / fs)
    return left, right, center


def qrs_correlation(sig, fs, meta0, meta1):

    l0, r0, c0 = get_qrs_bounds(meta0, fs)
    l1, r1, c1 = get_qrs_bounds(meta1, fs)
    leftw = min(c0-l0, c1-l1)
    rightw = min(r0-c0, r1-c1)

    qrs0 = sig[c0-leftw:c0+rightw, :]
    qrs1 = sig[c1-leftw:c1+rightw, :]

    cc = correlate(qrs0, qrs1)

    return cc


def classification(sig, hdr, metadata):

    qrs_classes = {}
    num_cyc = len(metadata)
    corrmat = np.zeros(num_cyc, num_cyc)
    fs = hdr["fs"]


    for i, meta0 in enumerate(metadata):

        for j in range(i+1, num_cyc):
            meta1 = metadata[j]

            corrmat[i, j] = qrs_correlation(sig, fs, meta0, meta1)







