#!/usr/bin/anv python
# coding: utf-8

import wfdb
import sys
import numpy as np

from argparse import ArgumentParser
from matplotlib import pyplot as plt
from sigbind import build_comb_filter, mean_spectrum, mains_filter, signal_channels
from qrsdetect import *


def build_args():
    parser = ArgumentParser()

    parser.add_argument(
        '-t', '--time-range',
        type=int,
        default=3,
        help='Time in seconds to display'
    )

    options, filenames = parser.parse_known_args()

    if not filenames:
        print("At least one input file should be specified")
        sys.exit(1)

    return options, filenames


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


def show_qrs(recordname, tend):

    # загрузили сигнал
    sig, fields = wfdb.rdsamp(recordname)
    sig = np.transpose(sig)
    fs = fields["fs"]
    num_chans = sig.shape[0]
    print(sig.shape)

    # следующие две строчки для наглядности. В реальной прорамме достаточно
    # вызвать
    # meta, foo = qrs_detection(sig, fs)

    ptstat = qrs_preprocessing(sig, fs)
    meta, strobe = qrs_detection(
        sig,
        fs=fs,
        bias=np.array([0.0]*num_chans),
        gain=np.array([1.0]*num_chans))

    # номера начального и конечного отсчета для отображения
    N1 = 0
    N2 = int(tend*fs)

    tt = np.linspace(float(N1)/fs, float(N2-1)/fs, N2-N1)

    plt.style.use("ggplot")
    plt.rcParams["figure.facecolor"] = "white"
    fig, axarr = plt.subplots(num_chans+1, 1, sharex=True)

    for chan, signal in enumerate(signal_channels(sig)):

        axarr[chan].plot(tt, signal[N1:N2])

        # пометим R-зубцы
        rx = []
        ry = []
        for x in meta:
            if tt[0] <= x["r_wave_center"] <= tt[-1]:
                rx.append(x["r_wave_center"])
                ry.append(signal[int(rx[-1]*fs)])

        axarr[chan].scatter(rx, ry, c="r")
        axarr[chan].plot(tt,strobe[N1:N2])

    # решающая статистика и стробы QRS-комплексов
    axarr[num_chans].plot(tt,ptstat[N1:N2])
    # axarr[num_chans].plot(tt,strobe[N1:N2])

    plt.show()


def main():
    options, filenames = build_args()

    #show_spectrums(sys.argv[1])
    show_qrs(filenames[0], options.time_range)


if __name__ == "__main__":
    main()