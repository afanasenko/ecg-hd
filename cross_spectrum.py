# coding: utf-8

"""
норма
s171012_132635
+ 0 (450 с) - 35/65 норма
+ 451 (290 с) - 17/132
+ 741 (300 с) - 70/33

"""

import sys

from argparse import ArgumentParser
from matplotlib import pyplot as plt

from zetlabreader import anaread
from sigbind import *
from sigsegment import extract_short_peaks


def build_args():
    parser = ArgumentParser()

    parser.add_argument(
        "-i", "--interval",
        type=int,
        default=150,
        help="Expected minimum R-peak interval"
    )

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
        "-w", "--window",
        type=int,
        default=0,
        help="Estimation window (seconds)"
    )

    parser.add_argument(
        "-o", "--out-file",
        type=str,
        help="Output file name"
    )

    return parser.parse_known_args()


def time_slice(signal, start_sec, end_sec, fs, ):
    start = np.round(int(fs * start_sec))
    end = np.round(int(fs * end_sec))

    if end:
        return signal[start:end]
    else:
        return signal[start:]


def main():
    options, files = build_args()
    if len(files) != 2:
        print("2 files expected")
        sys.exit(1)

    siglo, fs_lo, signal_name = anaread(files[0])
    print("Signal [{}], {} Hz loaded".format(len(siglo), fs_lo))

    sighi, fs_hi, signal_name = anaread(files[1])
    print("Signal [{}], {} Hz loaded".format(len(siglo), fs_hi))

    data = (
        time_slice(
            siglo - np.mean(siglo),
            options.offset_first,
            options.offset_first + options.window,
            fs_lo
        ),
        time_slice(
            sighi - np.mean(sighi),
            options.offset_first,
            options.offset_first + options.window,
            fs_hi
        )
    )

    plt.style.use("ggplot")
    plt.rcParams["font.family"] = "Verdana"
    f, axarr = plt.subplots(1, 1)

    peak_interval_ms = options.interval
    bias_ms = 1.9 * peak_interval_ms

    rpeaks, bks, hfs = extract_short_peaks(
        data[0],
        fs_lo,
        bias_window_ms=bias_ms,
        peak_interval_ms=peak_interval_ms
    )

    if len(rpeaks):
        print("{} cycles found".format(len(rpeaks)))
    else:
        print("No peaks")
        sys.exit(1)

    deltas = np.diff(rpeaks)
    q = np.percentile(deltas, [10, 50, 90])
    width = q[1]

    print("average cycle length: {}".format(width / fs_lo))

    rpk_hi = np.round(fs_hi * rpeaks / fs_lo).astype(int)

    sp = synchro_spectrum_r(
        data[1],
        rpk_hi,
        halfw=int(width/2),
        log_output=True
    )

    npoints = len(sp)
    xf = np.linspace(0.0, 0.5 * fs_hi, npoints)


    if options.out_file:
        with open(options.out_file, "w") as fo:
            for x, y in zip(xf, sp):
                fo.write("{}\t{}\n".format(x, y))
    else:
        axarr.plot(xf, sp)
        axarr.set_xlim([0, fs_hi/4])
        axarr.set_xlabel("Частота, Гц")
        axarr.set_ylabel("Амплитуда, дБ")
        axarr.grid(True)

        print("Look at the plots")
        plt.show()


if __name__ == "__main__":
    main()