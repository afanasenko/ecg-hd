# coding: utf-8

import sys
from matplotlib import pyplot as plt


def spload(filename):
    xf = []
    sp = []
    with open(filename, "r") as f1:
        for line in f1:
            p = [float(x.strip()) for x in line.split("\t")]
            xf.append(p[0])
            sp.append(p[1])

    return xf, sp

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(1)

    xf1, sp1 = spload(sys.argv[1])
    xf2, sp2 = spload(sys.argv[2])

    plt.style.use("ggplot")
    plt.rcParams["font.family"] = "Verdana"
    f, axarr = plt.subplots(1, 1)

    axarr.plot(xf1, sp1, '--', linewidth=2)
    axarr.plot(xf2, sp2, '-', linewidth=1)
    # axarr.set_xlim([0, fs_hi / 4])
    axarr.set_xlabel("Частота, Гц")
    axarr.set_ylabel("Амплитуда, дБ")
    axarr.grid(True)

    print("Look at the plots")
    plt.show()