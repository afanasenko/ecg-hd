# coding: utf-8

import numpy as np

from argparse import ArgumentParser
from sigsegment import extract_peaks_morpho, extract_short_peaks
from zetlabreader import anaread


def build_args():
    parser = ArgumentParser()

    parser.add_argument(
        "-i", "--interval",
        type=int,
        default=150,
        help="Expected minimum R-peak interval"
    )

    parser.add_argument(
        "-o", "--offset",
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

    return parser.parse_known_args()


def mean_heart_rate(signal, fs, peak_interval_ms, verbose=True):

    bias_ms = 1.9*peak_interval_ms

    rpeaks, bks, hfs = extract_short_peaks(signal, fs, bias_window_ms=bias_ms,
                                           peak_interval_ms=peak_interval_ms)

    if verbose:
        print("{} R-R intervals found".format(len(rpeaks)))

    if len(rpeaks):
        deltas = np.diff(rpeaks)
        q = np.percentile(deltas, [10, 50, 90])
        width = q[1]

        if verbose:
            print("Average R-R = {} samples".format(width))
            print("Confidence = [{:.1f}, {:.1f}]".format(
                60.0 * fs / q[2],
                60.0 * fs / q[0]
            ))

        return 60.0 * fs / width
    else:
        return 0


def main():
    options, files = build_args()

    for fn in files:
        data, fs, signal_name = anaread(fn)
        print("Signal [{}], {} Hz loaded".format(len(data), fs))

        start = np.round(int(fs*options.offset))
        end = np.round(int(fs*(options.offset + options.window)))

        if end:
            data = data[start:end]
        else:
            data = data[start:]

        beat_rate = mean_heart_rate(
            data,
            peak_interval_ms=options.interval,
            fs=fs
        )

        print("{}: {:.1f} bpm".format(
            fn,
            beat_rate
        ))


if __name__ == "__main__":
    main()