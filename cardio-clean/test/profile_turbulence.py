# coding: utf-8

import cProfile as profile
import pstats
import json

from cardio_clean.turbulence import turbulence_analyse

pr = profile.Profile()

metaname = "/Users/arseniy//Downloads/Test20191007.ecg.json"

print("Load...")
meta = json.load(open(metaname))
print("{} cycles".format(len(meta)))
print("Start...")

pr.enable()
turb_data, trend_data = turbulence_analyse(meta)
pr.disable()

print("TO={}, TS={}".format(trend_data["TO"], trend_data["TS"]))

pr.dump_stats('profile.pstat')

# Выводим 10 самых тяжелых функций
p = pstats.Stats('profile.pstat')
p.sort_stats('time').print_stats(10)

