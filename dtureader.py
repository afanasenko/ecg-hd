# coding: utf-8

import sys
import numpy as np
from matplotlib import pyplot as plt, rc


def dturead(filename):
    """
    Чтение записи, сделанной многоканальным осциллографом
    :param filename: полное имя файла
    :return: массив отсчетов (в физ. единицах), частота дискретизации, метка
    """

    if not filename.endswith(".dtu"):
        print("Расширение файла, отличное от dtu. Возможно неправильное "
              "чтение")

    headerlines = 8

    label = ""
    timestamps = []
    datablob = []

    with open(filename, "r") as fp:
        for i, line in enumerate(fp):
            if i == 0:
                label = line.strip().decode("utf-8")

            if i >= headerlines:
                try:
                    p = [float(x.strip()) for x in line.split()]
                except ValueError:
                    print("Ошибка в строке {} (нечисловые данные)\n{}".format(
                        i,
                        line
                    ))
                    break
                timestamps.append(p[0])

                datablob.append(p[1:])

    fs = 1.0 / (timestamps[1] - timestamps[0])

    d = np.array(datablob)
    print(d.shape)

    return d, fs, label


if __name__ == "__main__":

    if len(sys.argv) > 1:
        fn = sys.argv[1]
    else:
        fn = "Не задано имя файла"

    d, fs, n = dturead(fn)

    t = np.arange(0, len(d)/fs, 1.0/fs)

    rc('font', family="Verdana")
    plt.plot(t, d[:,0])
    #plt.xlim([540, 545])
    plt.title(n)

    print("Look at the plots")
    plt.show()