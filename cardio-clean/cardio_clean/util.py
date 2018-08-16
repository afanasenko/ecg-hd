# coding: utf-8

import wfdb
import numpy as np
from cardio_clean.cardioproc_api import read_buffer

def ecgread(filename):
    if filename.endswith(".ecg"):
        with open(filename, "rb") as fi:
            hdr, data = read_buffer(fi)
            return data, hdr

    else:
        data, fields = wfdb.rdsamp(filename)
        # rdsamp возвращает сигнал без смещения в физических единицах
        numch = data.shape[1]
        hdr = {
            "fs": fields["fs"],
            "adc_gain": np.array([1.0]*numch),
            "baseline": np.array([0.0]*numch),
            "samples": data.shape[0],
            "channels": data.shape[1]
        }

        return data, hdr