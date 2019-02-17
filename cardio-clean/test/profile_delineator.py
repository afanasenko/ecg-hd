# coding: utf-8

import cProfile as profile
import pstats

from cardio_clean.wavdetect import find_points
from cardio_clean.qrsdetect import qrs_detection
from cardio_clean.metadata import metadata_postprocessing
from cardio_clean.arrythmia import define_rythm
from cardio_clean.util import ecgread

pr = profile.Profile()

sig, header = ecgread("/Users/arseniy/SERDECH/data/ROXMINE/Rh2021")
fs = header["fs"]
if fs != 250:
    print("Warning! fs={}".format(fs))

print('start qrs_detection')
metadata = qrs_detection(sig,
                         fs=header["fs"],
                         minqrs_ms=20)[0]

print('start find_points')
find_points(sig,
            fs=header["fs"],
            metadata=metadata,
            bias=header["baseline"],
            gain=header["adc_gain"],
            debug=False)

print('start metadata_postprocessing')
metadata_postprocessing(metadata,
                        sig,
                        header)

print('start define_rythm')
pr.enable()
rithms = define_rythm(metadata)
pr.disable()

pr.dump_stats('profile.pstat')

# Выводим 10 самых тяжелых функций
p = pstats.Stats('profile.pstat')
p.sort_stats('time').print_stats(10)

