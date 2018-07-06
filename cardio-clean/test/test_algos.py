# coding: utf-8

from cardio_clean.metadata import *


def test_conv():
    fs = 250
    smp = fs

    assert ms_to_samples(samples_to_ms(smp, fs), fs) == smp
    assert samples_to_ms(fs, fs) == 1000.0
