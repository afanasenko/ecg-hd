# coding: utf-8

import numpy as np
from scipy.signal import lfilter

def dummy_shift(x, n):
    return np.concatenate((x[n:], np.ones(n)*x[-1]))


def qrs_preprocessing(sig, fs):
    """
    Предварительная обработка для выделения QRS
    :param sig: ЭКС (одноканальный или многоканальный)
    :param fs: частота дискретизации
    :return: характеристическая функция для дальнейшего детектирования
    """
    #TODO: реализовать синтез фильтров для произвольной fs
    if fs != 250:
        print("WARNING! выделение QRS для частоты дискретизации 250 Гц")

    result = None

    # чтобы одинаково обрабатывать одноканальгные и многоканальные сигналы,
    # добавляем размерность
    if len(sig.shape) == 1:
        sig = np.expand_dims(sig, axis=1)

    for channel in range(sig.shape[1]):
        # НЧ фильтр (1 - z ^ -6) ^ 2 / (1 - z ^ -1) ^ 2
        b = np.array([1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1], float)
        a = np.array([1, -2, 1], float)

        lp = lfilter(b, a, sig[:, channel])
        lp = dummy_shift(lp, 6)

        # слабый ВЧ фильтр z ^ -16 - [(1 - z ^ -32) / (1 - z ^ -1)]
        b = np.array([
            -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, -32,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
        ], float)
        a = np.array([1, -1], float)

        hp = lfilter(b, a, lp)
        hp = dummy_shift(hp, 16)

        # еще один ВЧ фильтр (производная)
        h = np.array([-1, -2, 0, 2, 1], float) / 8

        # сразу возводим в квадрат
        de = lfilter(h, [1.0], hp)**2
        de = dummy_shift(de, 2)

        # усреднение
        sm = lfilter(np.ones(31, float), [1.0], de)
        sm = dummy_shift(sm, 15)

        # решающие статистики для всех каналов просто суммируются
        if result is None:
            result = sm
        else:
            result += sm

    return sm / max(sm)