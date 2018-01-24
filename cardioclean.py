#!/usr/bin/env python
# coding: utf-8

import sys
import os
import yaml
import argparse
import wfdb
from sigbind import *


def build_options():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-c', '--config',
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'config.yaml'
        ),
        type=str,
        help='File name for associated yaml config'
    )
    parser.add_argument(
        '-i', '--input-file',
        required=True,
        type=str,
        help='Input file name (basename or full path to .dat / .hea file)'
    )

    parser.add_argument(
        '-o', '--output-file',
        required=True,
        type=str,
        help='Output file name (basename)'
    )

    options, other_params = parser.parse_known_args()
    return options


def main():
    options = build_options()

    try:
        with open(options.config, "r") as f:
            config = yaml.load(f)
    except Exception as e:
        print("Config file {} not found or invalid".format(options.config))
        sys.exit(1)

    recordname = ".".join(options.input_file.split(".")[:-1])

    sig, fields = wfdb.rdsamp(recordname)

    # число каналов берем из данных, а не из заголовка
    for channel in range(sig.shape[1]):

        if config["PREPROCESSING"].get("baseline_correction", True):
            ubsig = fix_baseline(sig[:, channel], fields["fs"], config[
                "BASELINE"]["unbias_window_ms"])
        else:
            ubsig = np.array(sig[:, channel], "float")

        if config["PREPORCESSING"].get("mains_correction", True):
            umsig = mains_filter(
                ubsig,
                fs=fields["fs"],
                mains=config["MAINS_FILTER"]["base_freq"],
                attenuation=config["MAINS_FILTER"]["attenuation"],
                aperture=config["MAINS_FILTER"]["fft_size"]
            )
        else:
            umsig = ubsig.copy()

        # восстановить исходный формат данных ??
        sig[:, channel] = umsig


if __name__ == "__main__":
    main()
