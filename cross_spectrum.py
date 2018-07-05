#!/usr/bin/env python
# coding: utf-8

import sys

from argparse import ArgumentParser
from matplotlib import pyplot as plt
from zetlabreader import anaread
import numpy as np

def build_args():
    parser = ArgumentParser()

    parser.add_argument(
        "-a", "--offset-first",
        type=int,
        default=0,
        help="Offset from start (seconds)"
    )

    parser.add_argument(
        "-b", "--offset-second",
        type=int,
        default=0,
        help="Offset from start (seconds)"
    )

    parser.add_argument(
        "-c", "--window",
        type=int,
        default=0,
        help="Estimation window (seconds)"
    )

    parser.add_argument(
        "-o", "--out-file",
        type=str,
        help="Output file name"
    )

    parser.add_argument(
        "-t", "--time-plot",
        action='store_true',
        default=False,
        help="Построение осциллограмм"
    )

    parser.add_argument(
        "-f", "--freq-plot",
        action='store_true',
        default=False,
        help="Построение спектров"
    )

    parser.add_argument(
        "-w", "--wavelet-plot",
        action='store_true',
        default=False,
        help="Построение вейвлет-диаграмм"
    )

    parser.add_argument(
        "-g", "--hist-plot",
        action='store_true',
        default=False,
        help="Построение гистограмм"
    )

    return parser.parse_known_args()


def crossplot_scope(sig1, fs1, sig1_name, sig2, fs2, sig2_name, out_name):

    f, axarr = plt.subplots(2, 1, sharex="True")

    tp1 = np.linspace(0, len(sig1)/fs1, fs1)
    axarr[0].plot(tp1, sig1)
    axarr[0].set_xlabel(u"Время, с")
    axarr[0].set_ylabel(sig1_name)
    axarr[0].grid(True)

    tp2 = np.linspace(0, len(sig2)/fs2, fs2)
    axarr[1].plot(tp2, sig2)
    axarr[1].set_xlabel(u"Время, с")
    axarr[1].set_ylabel(sig2_name)
    axarr[1].grid(True)

    if out_name:
        plt.savefig(out_name)
    else:
        plt.show(block=False)


def simple_spectrum(sig, fs):
    sp = mean_spectrum(
        sig,
        aperture=1024,  # окно БПФ
        log_output=True # Логарифмическая шкала амплитуд
    )

    npoints = len(sp)
    xf = np.linspace(0.0, 0.5 * fs, npoints) # шкала реальных частот
    return xf, sp


def crossplot_spectrum(sig1, fs1, sig1_name, sig2, fs2, sig2_name, out_name):

    xf1, sp1 = simple_spectrum(sig1, fs1)
    xf2, sp2 = simple_spectrum(sig2, fs2)

    plt.plot(xf1, sp1)
    plt.plot(xf2, sp2)
    plt.legend((sig1_name, sig2_name))
    plt.xlabel(u"Частота, Гц")
    plt.ylabel(u"Амплитуда, дБ")
    plt.grid(True)

    if out_name:
        plt.savefig(out_name)
    else:
        plt.show(block=False)


def time_slice(signal, start_sec, end_sec, fs):
    start = np.round(fs * start_sec).astype(int)

    if start >= len(signal):
        print("Начальное смещение превышает длительность сигнала")
        return []

    end = np.round(fs * end_sec).astype(int)

    if end < len(signal):
        print("Временное окно выходит за пределы длительности сигнала")
        return signal[start:end]
    else:
        return signal[start:]


def main():
    options, files = build_args()
    if len(files) != 2:
        print("2 files expected")
        sys.exit(1)

    sig1, fs1, sig1_name = anaread(files[0])
    print("1: запись {}, длительность [{}] с, fд={} Гц".format(
        sig1_name,
        len(sig1)/fs1,
        fs1
    ))

    sig2, fs2, sig2_name = anaread(files[1])
    print("1: запись {}, длительность [{}] с, fд={} Гц".format(
        sig2_name,
        len(sig2)/fs2,
        fs2
    ))

    sig1 = time_slice(
        sig1 - np.mean(sig1),
        options.offset_first,
        options.offset_first + options.window,
        fs1
    )

    sig2 = time_slice(
        sig2 - np.mean(sig2),
        options.offset_second,
        options.offset_second + options.window,
        fs2
    )

    outname = ""

    plt.style.use("ggplot")
    plt.rcParams["font.family"] = "Verdana"

    # Построение осциллограмм входных сигналов
    if options.time_plot:
        crossplot_scope(sig1, fs1, sig1_name, sig2, fs2, sig2_name, outname)

    # Построение спектров
    if options.freq_plot:
        crossplot_spectrum(sig1, fs1, sig1_name, sig2, fs2, sig2_name, outname)



if __name__ == "__main__":
    main()