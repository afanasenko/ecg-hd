#!/usr/bin/anv python
# coding: utf-8

import wfdb
import sys
import numpy as np

from matplotlib import pyplot as plt
from sigbind import build_comb_filter, mean_spectrum, mains_filter


def main():

    recordname = sys.argv[1]

    sig, fields = wfdb.rdsamp(recordname)
    sp = mean_spectrum(sig[:, 0])
    npoints = len(sp)
    xf = np.linspace(0.0, 0.5 * fields["fs"], npoints)

    fft_apert = 512
    mains_att = 0.05
    x, y = build_comb_filter(fields["fs"], fft_apert, mains_att)

    umsig = mains_filter(
        sig[:, 0],
        fs=fields["fs"],
        mains=50.0,
        attenuation=mains_att,
        aperture=fft_apert
    )

    sp2 = mean_spectrum(umsig)

    plt.rcParams["figure.facecolor"] = "white"

    fig, axarr = plt.subplots(2, 2)
    axarr[0,0].plot(xf, sp)
    axarr[1,0].plot(xf, sp2)

    happ = int(1 + fft_apert/2)
    axarr[0,1].plot(x[:happ], y[:happ], 'r')

    plt.show()


if __name__ == "__main__":
    main()