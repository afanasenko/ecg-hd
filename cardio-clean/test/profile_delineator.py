# coding: utf-8

import cProfile as profile

from cardio_clean.wavdetect import find_points
from cardio_clean.qrsdetect import qrs_detection
from demo_analysis import ecgread

pr = profile.Profile()

sig, header = ecgread("/Users/arseniy/SERDECH/data/ROXMINE/Rh1011")
fs = header["fs"]
if fs != 250:
    print("Warning! fs={}".format(fs))

s = sig[:, 0]

metadata, debugdata = qrs_detection(
    sig[:, :],
    fs=header["fs"],
    bias=header["baseline"],
    gain=header["adc_gain"],
    minqrs_ms=20)

pr.enable()
newmeta = find_points(s, header["fs"], metadata)
pr.disable()

pr.dump_stats('profile.pstat')
