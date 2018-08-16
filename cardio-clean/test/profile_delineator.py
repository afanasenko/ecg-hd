# coding: utf-8

import cProfile as profile
import pstats

from cardio_clean.wavdetect import find_points
from cardio_clean.qrsdetect import qrs_detection
from cardio_clean.util import ecgread

pr = profile.Profile()

sig, header = ecgread("/Users/arseniy/SERDECH/data/ROXMINE/Rh1011")
fs = header["fs"]
if fs != 250:
    print("Warning! fs={}".format(fs))

metadata = qrs_detection(
    sig,
    fs=header["fs"],
)[0]

pr.enable()
find_points(sig[:, 0], header["fs"], metadata, header["baseline"])
pr.disable()

pr.dump_stats('profile.pstat')

# Выводим 10 самых тяжелых функций
p = pstats.Stats('profile.pstat')
p.sort_stats('time').print_stats(10)

