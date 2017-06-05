#coding: utf-8

import sys
import os


import numpy as np
import matplotlib.pyplot as plt

root = os.path.dirname(os.path.abspath(__file__))
files = [os.path.join(root, x) for x in os.listdir(root) if x.endswith("features.tsv")]
print("\n".join(files))

plt.style.use("ggplot")
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size
fig, axarr = plt.subplots(1, 2)


ramp = []
tamp = []
rdur = []
tdur = []

spec = ['r', 'g', 'b', 'c', 'm', 'y', 'k']


for i, fn in enumerate(files):
    with open(fn, "r") as fo:
        ramp = []
        tamp = []
        rdur = []
        tdur = []
        for line in fo:
            p = [x.strip() for x in line.split("\t")]

            ramp.append(float(p[0]))
            tamp.append(float(p[2]))
            rdur.append(float(p[1]))
            tdur.append(float(p[3]))

        axarr[0].scatter(ramp, tamp, alpha=0.12, c=spec[i%len(spec)])
        axarr[0].set_xlabel("амплитуда R, мВ")
        axarr[0].set_ylabel("амплитуда T, мВ")

        axarr[1].scatter(tdur, tamp, alpha=0.1, c=spec[i%len(spec)])
        axarr[1].set_ylabel("амплитуда T, мВ")
        axarr[1].set_xlabel("длительность T, c")

plt.show()