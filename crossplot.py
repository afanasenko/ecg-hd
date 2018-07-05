#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np

from argparse import ArgumentParser
from matplotlib import pyplot as plt, rc
from zetlabreader import anaread
from sigbind import mean_spectrum, mains_filter
from scipy import signal


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

    parser.add_argument(
        "-m", "--mains-filter",
        action='store_true',
        default=False,
        help="Подавление сетевой помехи"
    )

    return parser.parse_known_args()


def crossplot_scope(sig1, fs1, sig1_name, sig2, fs2, sig2_name, out_name):
    if len(sig2):
        f, axarr = plt.subplots(2, 1, sharex="col")
    else:
        f, ax = plt.subplots(1, 1)
        axarr = [ax, None]

    tp1 = np.linspace(0, len(sig1) / fs1, len(sig1))
    axarr[0].plot(tp1, sig1)
    axarr[0].set_ylabel(sig1_name)
    axarr[0].grid(True)

    if len(sig2):
        tp2 = np.linspace(0, len(sig2) / fs2, len(sig2))
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
        log_output=True  # Логарифмическая шкала амплитуд
    )

    npoints = len(sp)
    xf = np.linspace(0.0, 0.5 * fs, npoints)  # шкала реальных частот
    return xf, sp


def crossplot_spectrum(sig1, fs1, sig1_name, sig2, fs2, sig2_name, out_name):
    f, ax = plt.subplots(1, 1)
    xf1, sp1 = simple_spectrum(sig1, fs1)
    plt.plot(xf1, sp1)

    if len(sig2):
        xf2, sp2 = simple_spectrum(sig2, fs2)
        ax.plot(xf2, sp2)
        ax.legend((sig1_name, sig2_name))

    ax.set_xlabel(u"Частота, Гц")
    ax.set_ylabel(u"Амплитуда, дБ")
    ax.grid(True)

    if out_name:
        plt.savefig(out_name)
    else:
        plt.show(block=False)


def simple_histogram(sig):
    f, x = np.histogram(sig,bins=50)

    return x, np.array([f[0]] + list(f))


def crossplot_hist(sig1, fs1, sig1_name, sig2, fs2, sig2_name, out_name):

    f, ax = plt.subplots(1, 1)
    x1, f1 = simple_histogram(sig1)
    plt.step(x1, f1)

    if len(sig2):
        x2, f2 = simple_histogram(sig2)
        ax.step(x2, f2)
        ax.legend((sig1_name, sig2_name))

    ax.set_xlabel(u"Значение")
    ax.set_ylabel(u"Частота")
    ax.grid(True)

    if out_name:
        plt.savefig(out_name)
    else:
        plt.show(block=False)


def simple_wlt(sig, fs, signame, scales, ax):
    cwtmatr = signal.cwt(sig, signal.ricker, scales)

    ax.imshow(
        cwtmatr,
        extent=[0, int(len(sig) / fs), min(scales), max(scales)],
        cmap='nipy_spectral',
        aspect='auto',
        vmax=abs(cwtmatr).max(),
        vmin=-abs(cwtmatr).max()
    )

    ax.get_yaxis().set_visible(False)
    ax.set_xlabel(u"Время, с")
    ax.set_title(signame)


def crossplot_cwt(sig1, fs1, sig1_name, sig2, fs2, sig2_name, out_name):
    if len(sig2):
        f, axarr = plt.subplots(2, 1, sharex="col")
    else:
        f, ax = plt.subplots(1, 1)
        axarr = [ax, None]

    widths = [1]
    while widths[-1] < 500:
        widths.append(1.05 * widths[-1])

    simple_wlt(sig1, fs1, sig1_name, widths, axarr[0])

    if len(sig2):
        simple_wlt(sig2, fs2, sig2_name, widths, axarr[1])

    if out_name:
        plt.savefig(out_name)
    else:
        plt.show(block=False)


def time_slice(sig, start_sec, end_sec, fs):
    start = np.round(fs * start_sec).astype(int)

    if start >= len(sig):
        print("Начальное смещение превышает длительность сигнала")
        return []

    end = np.round(fs * end_sec).astype(int)

    if end < len(sig):
        return sig[start:end]
    else:
        print("Временное окно выходит за пределы длительности сигнала")
        return sig[start:]


def processing(file1, file2, use_mains_filter, show_scope, show_freq,
               show_hist, show_cwt, file_output, offset_first,
               offset_second, window):
    sig1, fs1, sig1_name = anaread(file1)
    if not sig1_name:  # если нет в описании
        sig1_name = u"первый"
    print(u"1: запись \"{}\", длительность {} с, fд={} Гц".format(
        sig1_name,
        len(sig1) / fs1,
        fs1
    ))

    if file2:
        sig2, fs2, sig2_name = anaread(file2)
        if not sig2_name:  # если нет в описании
            sig2_name = u"второй"
        print(u"2: запись \"{}\", длительность {} с, fд={} Гц".format(
            sig2_name,
            len(sig2) / fs2,
            fs2
        ))
    else:
        sig2 = []
        fs2 = 0
        sig2_name = u"не задан"
        print(u"Второй файл не задан")

    if use_mains_filter:
        sig1 = mains_filter(sig1, fs1, bias=np.mean(sig1), mains=50.0,
                            attenuation=0.01, aperture=1024)

    sig1 = time_slice(
        sig1 - np.mean(sig1),
        offset_first,
        offset_first + window,
        fs1
    )

    if len(sig2):
        if use_mains_filter:
            sig2 = mains_filter(sig2, fs2, bias=np.mean(sig2), mains=50.0,
                            attenuation=0.01, aperture=1024)

        sig2 = time_slice(
            sig2 - np.mean(sig2),
            offset_second,
            offset_second + window,
            fs2
        )

    if not len(sig1) and not len(sig2):
        print("Невозможно выделить заданные временные интервалы")
        return

    plt.style.use("ggplot")
    plt.rcParams["font.family"] = "Verdana"
    rc('font', family="Verdana")

    # Построение осциллограмм входных сигналов
    if show_scope:
        if file_output:
            outname = file_output + u"_осциллограммы.png"
        else:
            outname = ""
        crossplot_scope(sig1, fs1, sig1_name, sig2, fs2, sig2_name, outname)

    # Построение спектров
    if show_freq:
        if file_output:
            outname = file_output + u"_спектры.png"
        else:
            outname = ""
        crossplot_spectrum(sig1, fs1, sig1_name, sig2, fs2, sig2_name, outname)

    # Построение гистограмм
    if show_hist:
        if file_output:
            outname = file_output + u"_гистограммы.png"
        else:
            outname = ""
        crossplot_hist(sig1, fs1, sig1_name, sig2, fs2, sig2_name, outname)

    # Построение вейвлетов
    if show_cwt:
        if file_output:
            outname = file_output + u"_вейвлеты.png"
        else:
            outname = ""
        crossplot_cwt(sig1, fs1, sig1_name, sig2, fs2, sig2_name, outname)


def main():
    options, files = build_args()
    if len(files) < 2:
        print("Укажите 2 файла для сравнения")
        sys.exit(1)

    processing(
        file1=files[0],
        file2=files[1],
        use_mains_filter=options.mains_filter,
        show_scope=options.time_plot,
        show_freq=options.freq_plot,
        show_hist=options.hist_plot,
        show_cwt=options.wavelet_plot,
        file_output=options.out_file,
        offset_first=options.offset_first,
        offset_second=options.offset_second,
        window=options.window,
    )

    # чтобы не позакрывались графики
    if not options.file_output:
        print("Нажмите любую клавишу для завершения")
        sys.stdin.read(1)



if __name__ == "__main__":
    main()