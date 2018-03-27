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


def qrs_detection(sig, fs, minqrs_ms=20):
    """

    :param sig: ЭКС (одноканальный или многоканальный)
    :param fs: частота дискретизации
    :param minqrs_ms: минимальная длительность QRS-комплекса
    :param debug: вывод дополнительных данных
    :return: qrs_metadata (список найденных комплексов), отладочные данные
    """

    pp = qrs_preprocessing(sig, fs)

    # 4-секундное окно для адаптивного порога
    thresh_wnd = 4*fs
    halfw = int(thresh_wnd/2)

    qrsmask = np.zeros(len(pp), int)

    inside = False
    qrs_start = 0
    qrs_num = 0
    minqrs_smp = fs * float(minqrs_ms) / 1000
    pilot_channel = 0

    qrs_metadata = []

    for i, v in enumerate(pp):
        i0 = max(0, i-halfw)
        i1 = min(len(pp), i+halfw)
        thresh = np.mean(pp[i0:i1])
        flag = 1 if v > thresh else 0
        qrsmask[i] = flag

        if flag and not inside:
            qrs_start = i
            inside = True

        if not flag and inside:
            qrs_end = i
            inside = False
            qrs_len = qrs_end - qrs_start
            if qrs_len >= minqrs_smp:

                #TODO: R-зубец искать в одном канале или во всех?
                rpk = qrs_start + np.argmax(
                    sig[qrs_start:qrs_end, pilot_channel])

                qrs_metadata.append({
                    "cycle_num": qrs_num,
                    "qrs_start": float(qrs_start)/fs,
                    "qrs_end": float(qrs_end)/fs,
                    "q_wave_center": None,
                    "r_wave_center": float(rpk)/fs,
                    "s_wave_center": None
                })
                qrs_num += 1

    return qrs_metadata, qrsmask
