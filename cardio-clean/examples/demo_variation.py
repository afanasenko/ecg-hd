# coding: utf-8


import os
from matplotlib import pyplot as plt
import json


from cardio_clean.metadata import metadata_postprocessing, calculate_histogram
from cardio_clean.wavdetect import find_points
from cardio_clean.spectralvariation import *
from cardio_clean.statvariation import rhythm_stats

from cardio_clean.sigbind import fix_baseline
from cardio_clean.qrsdetect import qrs_detection
from cardio_clean.util import ecgread, signal_channels


def show_rspec(filename, chan, smp_from=0, smp_to=0):
    sig, header = ecgread(filename)

    print(sig.shape)

    fs = header["fs"]

    sig = fix_baseline(
        sig,
        fs=fs,
        bias_window_ms=1500
    )

    if smp_to:
        smp_to = min(smp_to, sig.shape[0])
    else:
        smp_to = sig.shape[0]

    metadata, foo = qrs_detection(
        sig[smp_from:smp_to,:],
        fs=header["fs"]
    )

    find_points(
        sig[smp_from:smp_to, :],
        fs=header["fs"],
        bias=header["baseline"],
        gain=header["adc_gain"],
        metadata=metadata
    )

    metadata_postprocessing(
        metadata,
        sig[smp_from:smp_to, :],
        header
    )

    compl = {}
    for x in metadata:
        c = x["complex_type"]
        compl[c] = compl.get(c,0) + 1
    print(json.dumps(compl, indent=1))

    print("stat started")
    stat_vals = rhythm_stats(metadata)
    print(json.dumps(stat_vals, indent=1))

    print("spectrun started")
    v, r, fp, sp = rhythm_spectrum(metadata)
    print("df = {:.5f} Hz".format(fp[1]))
    print(v)
    plt.plot(fp, sp, "b")

    print("Look at the plots")
    plt.show()


def main():
    # Rh2022 = qr, noise
    # Rh2021 - Rs, extracyc
    # Rh2024 - p(q)RsT
    # Rh2025 = rs
    # Rh2010 - дрейф, шум, артефакты
    # 2004 av block

    #filename = "/Users/arseniy/SERDECH/data/PHYSIONET/I11"
    #filename = "/Users/arseniy/SERDECH/data/PHYSIONET/104"
    filename = "/Users/arseniy/SERDECH/data/PHYSIONET/222"
    #filename = "../test/testI59.ecg"
    #filename = "../test/TestFromDcm.ecg"
    #filename = "/Users/arseniy/SERDECH/data/ROXMINE/Rh2004"

    if not filename.endswith(".ecg") and not os.path.isfile(filename + ".hea"):
        print("Файл не найден")
        return

    show_rspec(
        filename,
        chan=1,  # common_signal_names.index("I")
    )


if __name__ == "__main__":
    main()