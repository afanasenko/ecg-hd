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

from sigbind import fix_baseline, mains_filter


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


def logsetup():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)s [%(name)s.%(module)s.%(funcName)s:%(lineno)d] %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S")
    # log to a stream
    stream_handler = logging.StreamHandler(stream=sys.stdout)
    # stream_handler.setLevel(logging.INFO)
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(formatter)
    # log to file
    log_filename = '{}.log'.format(
        os.path.splitext(os.path.basename(__file__))[0])
    file_handler = logging.FileHandler(log_filename, encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    # add handlers
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    # capture warnings from other modules that don't use logging, but use module 'warnings'
    logging.captureWarnings(capture=True)
    warnings_logger = logging.getLogger('py.warnings')
    warnings_logger.addHandler(file_handler)
    warnings_logger.addHandler(stream_handler)


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

    sig, fields = wfdb.rdsamp(src_file)

    # число каналов берем из данных, а не из заголовка
    for channel in range(sig.shape[1]):

        if config["PREPROCESSING"].get("baseline_correction", True):
            ubsig = fix_baseline(sig[:, channel], fields["fs"], config[
                "BASELINE"]["unbias_window_ms"])
        else:
            ubsig = np.array(sig[:, channel], "float")

        if config["PREPROCESSING"].get("mains_correction", True):
            umsig = mains_filter(
                ubsig,
                fs=fields["fs"],
                mains=config["MAINS_FILTER"]["base_freq"],
                attenuation=config["MAINS_FILTER"]["attenuation"],
                aperture=config["MAINS_FILTER"]["fft_size"]
            )
        else:
            umsig = ubsig.copy()

        sig[:, channel] = umsig

    wfdb.wrsamp(
        dest_file,
        fs=fields["fs"],
        units=fields["units"],
        sig_name=fields["sig_name"],
        comments=fields["comments"],
        p_signal=sig,
        fmt=fields.get("fmt", ["16"]*sig.shape[1])
    )


if __name__ == "__main__":
    main()
