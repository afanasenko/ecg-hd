#!/usr/bin/anv python
# coding: utf-8

from matplotlib import pyplot as plt
from sigbind import *

def main():

    recordname = ""

    # sig, fields = wfdb.srdsamp(recordname)
    x, y = sigbind.build_comb_filter(250, 1024, 0.05)

    plt.rcParams["figure.facecolor"] = "white"

    fig, axarr = plt.subplots(1, 1)
    axarr.plot(x[:513], y[:513], 'r')

    plt.show()


if __name__ == "__main__":
    main()