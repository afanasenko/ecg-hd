#!/usr/bin/anv python
# coding: utf-8

import wfdb
import sys
import numpy as np

from matplotlib import pyplot as plt
from sigbind import build_comb_filter, mean_spectrum, mains_filter
from qrsdetect import *


def show_spectrums(recordname):

    sig, fields = wfdb.rdsamp(recordname)
    sp = mean_spectrum(sig[:, 0])
    npoints = len(sp)
    xf = np.linspace(0.0, 0.5 * fields["fs"], npoints)

    fft_apert = 512
    mains_att = 0.05
    x, y = build_comb_filter(fields["fs"], fft_apert, mains_att)

    umsig = mains_filter(
        sig[:, 0],
        fs=fields["fs"],
        mains=50.0,
        attenuation=mains_att,
        aperture=fft_apert
    )

    sp2 = mean_spectrum(umsig)

    plt.rcParams["figure.facecolor"] = "white"

    fig, axarr = plt.subplots(2, 2)
    axarr[0,0].plot(xf, sp)
    axarr[1,0].plot(xf, sp2)

    happ = int(1 + fft_apert/2)
    axarr[0,1].plot(x[:happ], y[:happ], 'r')

    plt.show()


def show_qrs(recordname):

    # загрузили сигнал
    sig, fields = wfdb.rdsamp(recordname)
    fs = fields["fs"]

    # следующие две строчки для наглядности. В реальной прорамме достаточно
    # вызвать
    # meta, foo = qrs_detection(sig, fs)

    ptstat = qrs_preprocessing(sig, fs)
    meta, strobe = qrs_detection(sig, fs)

    # номера начального и конечного отсчета для отображения
    N1 = 0
    N2 = 900

    # номер канала
    chan = 0

    tt = np.linspace(float(N1)/fs, float(N2-1)/fs, N2-N1)

    plt.rcParams["figure.facecolor"] = "white"
    fig, axarr = plt.subplots(2, 1)

    axarr[0].plot(tt, sig[N1:N2, chan])

    # пометим R-зубцы
    rx = []
    ry = []
    for x in meta:
        if x["R-peak"] <= tt[-1]:
            rx.append(x["R-peak"])
            ry.append(sig[int(rx[-1]*fs), chan])

    axarr[0].scatter(rx, ry)

    # решающая статистика и стробы QRS-комплексов
    axarr[1].plot(tt,ptstat[N1:N2])
    axarr[1].plot(tt,strobe[N1:N2])

    plt.show()


def main():
    #show_spectrums(sys.argv[1])
    show_qrs("/Users/arseniy/Downloads/ROXMINE/Rh2010")


if __name__ == "__main__":
    main()