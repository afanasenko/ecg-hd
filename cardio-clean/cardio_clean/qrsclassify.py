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


def get_qrs_bounds(meta, fs):
    left = int(meta["qrs_start"] * fs)
    right = int(meta["qrs_end"] * fs)
    center = int(meta["r_wave_center"][0] * fs)
    return left, right, center


def qrs_correlation(sig, fs, meta0, meta1):

    l0, r0, c0 = get_qrs_bounds(meta0, fs)
    l1, r1, c1 = get_qrs_bounds(meta1, fs)
    leftw = min(c0-l0, c1-l1)
    rightw = min(r0-c0, r1-c1)

    qrs0 = sig[c0-leftw:c0+rightw, :]
    qrs1 = sig[c1-leftw:c1+rightw, :]

    cc = correlate(qrs0, qrs1)

    return float(cc)


def extract_qrs(sig, fs, meta):
    left, right, center = get_qrs_bounds(meta, fs)
    return sig[left:right, :], center - left


def qrs_reference_correlation(ref, sig, offset):

    smp = ref.shape[0]

    cc = correlate(
        ref,
        sig[offset:offset+smp, :]
    )

    return cc


def classification(sig, hdr, metadata):

    num_cyc = len(metadata)
    corrmat = np.ones((num_cyc, num_cyc), float)
    fs = hdr["fs"]

    for i in range(num_cyc):
        for j in range(i+1, num_cyc):
            cc = qrs_correlation(sig, fs, metadata[i], metadata[j])
            corrmat[i,j] = cc

    print(np.min(corrmat))


def incremental_classifier(sig, hdr, metadata, classgen_t=0.9):
    """

    :param sig:
    :param hdr:
    :param metadata:
    :param classgen_t: нижний порог коэффициента корреляции на создание нового класса
    :return:
    """

    num_cyc = len(metadata)
    fs = hdr["fs"]

    first_qrs, first_c = extract_qrs(sig, fs, metadata[1])

    # инициализировать первый класс вторым комплексом, т.к. 1-й может быть
    # обрезан
    qrs_classes = [
        {
            "accumulator": first_qrs,
            "center": first_c,
            "samples": {1}
        }
    ]

    for i in range(num_cyc):

        l1, r1, c1 = get_qrs_bounds(metadata[i], fs)
        cormat = np.zeros(len(qrs_classes), float)
        for c, qrsc in enumerate(qrs_classes):

            cormat[c] = qrs_reference_correlation(
                qrsc["accumulator"] / len(qrsc["samples"]),
                sig,
                c1 - qrsc["center"]
            )

        max_cc = np.max(cormat)

        if max_cc < classgen_t:
            # создание нового класса
            new_qrs, new_c = extract_qrs(sig, fs, metadata[i])
            qrs_classes.append(
                {
                    "accumulator": new_qrs,
                    "center": new_c,
                    "samples": {i}
                }
            )
        else:
            # выбор существуюущего класса
            classid = np.argmax(cormat)
            qcl = qrs_classes[classid]

            # обновление усредненного цикла в выбранном классе
            left_acc = qcl["center"]
            right_acc = qcl["accumulator"].shape[0] - left_acc

            if c1 - left_acc >= 0 and c1 + right_acc <= sig.shape[0]:

                qcl["accumulator"] += sig[c1-left_acc:c1+right_acc, :]
                qcl["samples"].add(i)

    return qrs_classes







