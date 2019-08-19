# coding: utf-8

import cProfile as profile
import pstats

from cardio_clean.wavdetect import find_points
from cardio_clean.qrsdetect import qrs_detection
from cardio_clean.qrsclassify import incremental_classifier
from cardio_clean.metadata import metadata_postprocessing
from cardio_clean.arrythmia import define_rythm
from cardio_clean.util import ecgread

pr = profile.Profile()

sig, header = ecgread("/Users/arseniy/SERDECH/data/Holter_24h")

fs = header["fs"]
if fs != 250:
    print("Warning! fs={}".format(fs))

ns = 18000*fs

sig = sig[:ns,:]

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
pr.enable()
metadata_postprocessing(metadata,
                        sig,
                        header)
pr.disable()

print('start incremental_classifier')

qrs_classes = incremental_classifier(
    sig,
    header,
    metadata,
    classgen_t=0.7,
    include_data=3
)


print('start define_rythm')

rithms = define_rythm(metadata)

pr.dump_stats('profile.pstat')

# Выводим 10 самых тяжелых функций
p = pstats.Stats('profile.pstat')
p.sort_stats('time').print_stats(10)

