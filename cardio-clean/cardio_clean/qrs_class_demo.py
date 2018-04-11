#!/usr/bin/anv python
# coding: utf-8

import wfdb
import sys
import numpy as np

from argparse import ArgumentParser
from matplotlib import pyplot as plt
from sigbind import fix_baseline
from qrsdetect import *
from qrsclassify import *
from graphresults import ecgread


def build_args():
    parser = ArgumentParser()

    parser.add_argument(
        '-t', '--time-range',
        type=int,
        default=6,
        help='Time in seconds to display'
    )

    options, filenames = parser.parse_known_args()
    if not filenames:
        filenames.append("/Users/arseniy/SERDECH/data/ROXMINE/Rh2002")
        #filenames.append("/Users/arseniy/heart-research/cardio-clean/test
        # /TestFromDcm.ecg")


    if not filenames:
        print("At least one input file should be specified")
        sys.exit(1)

    return options, filenames


def get_qrsclass(recordname, tend):

    # загрузили сигнал
    sig, hdr = ecgread(recordname)
    fs = hdr["fs"]

    sig = fix_baseline(
        sig,
        fs=fs,
        bias_window_ms=1500
    )

    meta, foo = qrs_detection(
        sig,
        fs=fs,
        bias=hdr["baseline"],
        gain=hdr["adc_gain"]
    )

    qrs_classes = incremental_classifier(
        sig,
        hdr,
        meta,
        classgen_t=0.7
    )

    print("{} cycles found".format(len(meta)))
    print("{} classes detected".format(len(qrs_classes)))
    print("{} complexes not classified".format(len(
        [1 for x in meta if x.get("qrs_class_id", None) is None]
    )))

    show_classes = min(len(qrs_classes), 5)
    fig, axarr = plt.subplots(hdr["channels"], show_classes, sharex=True)

    for i, qcl in enumerate(qrs_classes):
        if i >= show_classes:
            break

        samples = qcl["average"].shape[0]
        t = np.linspace(0.0, 1000.0*samples/fs, samples)

        for chan, wave in signal_channels(qcl["average"]):

            axarr[chan, i].plot(t, wave)

            # подписи отведений
            if i == 0:
                axarr[chan,i].set_ylabel("I"*(chan+1), rotation=0, fontsize=16)

        axarr[0,i].set_title("[{}]".format(qcl["count"]))
        axarr[-1,i].set_xlabel("t, ms")

    plt.show()


def main():
    options, filenames = build_args()

    get_qrsclass(filenames[0], options.time_range)


if __name__ == "__main__":
    main()