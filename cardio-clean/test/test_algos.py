# coding: utf-8

from cardio_clean.metadata import *
from cardio_clean.pvcdetect import *


def test_conv():
    fs = 250
    smp = fs

    assert ms_to_samples(samples_to_ms(smp, fs), fs) == smp
    assert samples_to_ms(fs, fs) == 1000.0


def test_pvc_epi():

    test_data = [
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 1.0},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 1.0},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 0.45},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 0.45},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 0.75},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 0.44},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 0.47},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 0.77},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 1.0},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 1.0},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 1.0},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 1.0},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 1.0},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 1.0},
        {"complex_type": "V", "r_pos": [0], "flags": "", "RR": 1.0},
    ]

    fs = 250
    p = 0
    for x in test_data:
        x["r_pos"][0] = p
        p += x["RR"]*fs

    marks = detect_pvc_episodes(test_data, fs=fs)
    assert len([x for x in marks if x != ""]) == 4


if __name__ == "__main__":
    test_pvc_epi()