#!/usr/bin/anv python
# coding: utf-8

import wfdb
import sys
import numpy as np

from argparse import ArgumentParser
from matplotlib import pyplot as plt
from cardio_clean.sigbind import fix_baseline
from cardio_clean.qrsdetect import *
from cardio_clean.qrsclassify import *
from cardio_clean.wavdetect import find_points
from cardio_clean.util import ecgread


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
        #filenames.append("/Users/arseniy/SERDECH/data/ROXMINE/Rh2004")
        filenames.append("/Users/arseniy/SERDECH/data/PHYSIONET/I16")
        #filenames.append("/Users/arseniy/heart-research/cardio-clean/test
        # /TestFromDcm.ecg")
        # 1003 2018

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

    metadata = qrs_detection(
        sig,
        fs=fs
    )[0]

    find_points(
        sig[:, 0],
        fs=fs,
        metadata=metadata,
        bias=hdr["baseline"],
        gain=hdr["adc_gain"]
    )

    qrs_classes = incremental_classifier(
        sig,
        hdr,
        metadata,
        classgen_t=0.7,
        include_data=3
    )

    print("{} cycles found".format(len(metadata)))
    print("{} classes detected".format(len(qrs_classes)))
    print("{} complexes not classified".format(len(
        [1 for x in metadata if x.get("qrs_class_id", None) is None]
    )))

    show_classes = min(len(qrs_classes), 5)
    fig, axarr = plt.subplots(hdr["channels"], show_classes, sharex="col")

    for i, qcl in enumerate(qrs_classes):
        # ограничение на число отображаемых классов, если их слишком много
        if i >= show_classes:
            break

        # расчет временных точек
        samples = qcl["average"].shape[0]
        t = np.linspace(0.0, 1000.0*samples/fs, samples)

        # построение усредненного комплекса
        for chan, wave in signal_channels(qcl["average"]):

            if show_classes == 1:
                axarr[chan].plot(t, wave)
                axarr[chan].grid(True)
            else:
                axarr[chan, i].plot(t, wave)
                axarr[chan, i].grid(True)

            # подписи отведений
            if i == 0:
                if show_classes == 1:
                    axarr[chan].set_ylabel("I"*(chan+1), rotation=0,
                                             fontsize=16)
                else:
                    axarr[chan, i].set_ylabel("I" * (chan + 1), rotation=0,
                                              fontsize=16)

        if show_classes == 1:
            axarr[0].set_title("[{}]".format(qcl["count"]))
            axarr[-1].set_xlabel("t, ms")
        else:
            axarr[0,i].set_title("[{}]".format(qcl["count"]))
            axarr[-1,i].set_xlabel("t, ms")

    plt.show()


def main():
    options, filenames = build_args()

    get_qrsclass(filenames[0], options.time_range)


if __name__ == "__main__":
    main()