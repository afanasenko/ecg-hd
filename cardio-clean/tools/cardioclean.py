#!/usr/bin/env python
# coding: utf-8

import sys
import os
import yaml
import argparse
import wfdb
import logging
import shutil

import numpy as np

from cardio_clean.sigbind import fix_baseline, mains_filter


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

    if options.input_file == "stdin":
        for line in sys.stdin:
            src_file = line.strip()
            recordname = ".".join(src_file.split(".")[:-1])
            dest_file = os.path.basename(recordname) + \
                        config["BATCH_PROCESSING"]["output_file_suffix"]

            print("processing {}".format(recordname))
            process_one_record(recordname, dest_file, config)

            attr_file_in = os.path.join(recordname + ".atr")
            attr_file_out = dest_file + ".atr"

            if os.path.isfile(attr_file_in):
                shutil.copyfile(attr_file_in, attr_file_out)
    else:
        process_one_record(
            src_file=".".join(options.input_file.split(".")[:-1]),
            dest_file=options.output_file,
            config=config
        )


def process_one_record(src_file, dest_file, config):

    data, fields = wfdb.rdsamp(src_file)
    numch = data.shape[1]
    hdr = {
        "fs": fields["fs"],
        "adc_gain": np.array([1.0] * numch),
        "baseline": np.array([0.0] * numch),
        "samples": data.shape[0],
        "channels": data.shape[1]
    }

    if config["PREPROCESSING"].get("baseline_correction", True):
        ubsig = fix_baseline(
            data,
            fields["fs"],
            config["BASELINE"]["unbias_window_ms"]
        )
    else:
        ubsig = np.array(data, "float")

    if config["PREPROCESSING"].get("mains_correction", True):
        umsig = mains_filter(
            ubsig,
            fs=fields["fs"],
            bias=hdr["baseline"],
            mains=config["MAINS_FILTER"]["base_freq"],
            attenuation=config["MAINS_FILTER"]["attenuation"],
            aperture=config["MAINS_FILTER"]["fft_size"]
        )
    else:
        umsig = ubsig.copy()

    wfdb.wrsamp(
        dest_file,
        fs=fields["fs"],
        units=fields["units"],
        sig_name=fields["sig_name"],
        comments=fields["comments"],
        p_signal=umsig,
        fmt=fields.get("fmt", ["16"]*umsig.shape[1])
    )


if __name__ == "__main__":
    main()
