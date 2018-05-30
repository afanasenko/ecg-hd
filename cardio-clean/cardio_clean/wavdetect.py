# coding: utf-8

import numpy as np
from scipy.signal import convolve

from sigbind import signal_channels


def fexpand(x):
    """
    Прореживание вдвое
    :param x:
    :return:
    """
    n = len(x)
    y = np.array([0]*2*n, float)
    for i in range(len(y)):
        if i%2 == 0:
            y[i] = x[i/2]

    return y


def ddwt(x, num_scales=5):
    """
    Дискретное вейвлет-преобразование без прореживания
    :param x: входной сигнал
    :param num_scales: число уровней разложения
    :return:
    """
    h = np.array([1, 3, 3, 1], float) / 8
    g = np.array([2, -2], float)

    decomposition = []
    ap = x.copy()
    decomposition.append(ap.copy())  # на нулевом уровне храним исходный сигнал

    for s in range(num_scales):

        # TODO: сразу компенсировать задержку
        hif = convolve(ap, g, mode="same")
        decomposition.append(hif)
        ap = convolve(ap, h, mode="same")

        # вместо прореживания сигналов (Маллат) расширяем характеристики
        # фильтров
        h = fexpand(h)
        g = fexpand(g)

    return decomposition


def find_points(x):

    scales = ddwt(x)



