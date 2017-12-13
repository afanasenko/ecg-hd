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
    if len(sys.argv) < 5:
        sys.exit(1)

    xf1, sp1 = spload(sys.argv[1])

    tit1 = sys.argv[2]

    xf2, sp2 = spload(sys.argv[3])

    tit2 = sys.argv[4]

    f1 = max(int(sys.argv[5]), xf1[0])
    f2 = min(int(sys.argv[6]), xf1[-1])

    plt.style.use("ggplot")
    plt.rcParams["font.family"] = "Verdana"
    f, axarr = plt.subplots(1, 1)


    axarr.plot(xf1, sp1, '-', linewidth=1)
    axarr.plot(xf2, sp2, '--', linewidth=2)
    axarr.set_xlim([f1, f2])
    axarr.set_xlabel("Частота, Гц")
    axarr.set_ylabel("Амплитуда, дБ")
    axarr.legend([tit1, tit2])
    axarr.grid(True)

    print("Look at the plots")
    plt.show()