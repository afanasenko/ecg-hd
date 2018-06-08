# coding: utf-8

import json
import numpy as np
from scipy.signal import convolve, argrelmax, argrelmin


def zcfind(x, single=True, lb=0, rb=0):

    zc = []

    if rb:
        rb = min(rb, len(x))
    else:
        rb = len(x)

    for i in range(lb+1, rb):
        if x[i-1]*x[i] < 0:
            #w1 = float(abs(x[i-1]))
            #w2 = float(abs(x[i]))
            if single:
                return i-1

            zc.append(i-1)

    if single:
        return zc[0]
    else:
        return np.array(zc)


def wave_bounds(x, wav, limit):

    w1 = wav[0]
    w2 = wav[1]

    for i in range(wav[1], wav[1] + limit):
        if x[i-1]*x[i] < 0:
            w2 = i
            break

    for i in np.arange(limit):
        if x[w1-i] * x[w1-i-1] < 0:
            w1 -= i
            break

    return w1, w2


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


def pksearch(modes, derivative):

    # Удаляем ступеньки
    mbuf = []
    for pos, val in modes:
        if mbuf:
            # повторные с одним знаком
            if val * mbuf[-1][1] > 0:
                # оставляем максимальный
                if abs(val) > abs(mbuf[-1][1]):
                    mbuf[-1] = (pos, val)
                continue

        mbuf.append((pos, val))

    # кодируем найденные фронты знаками +/-
    signcode = "".join(["+" if x[1] > 0 else "-" for x in mbuf])

    params = {
        "qrsType": None,
        "qWavePositon": None,
        "rWavePositon": None,
        "sWavePositon": None,
        "qWaveHeight": 0,
        "rWaveHeight": 0,
        "sWaveHeight": 0,
    }

    if signcode == "-+":
        params["qrsType"] = "qs"


    # Фронты самого мощного (R или  qs) зубца
    maxpair = (0,0)
    for i, posval in enumerate(modes):
        if i:
            diff = abs(posval[1]) + abs(modes[i - 1][1])
            if diff > maxpair[1]:
                maxpair = (i-1, diff)

    i0 = maxpair[0]

    params["rWavePositon"] = zcfind(
        derivative,
        single=True,
        lb=modes[i0][0],
        rb=modes[i0 + 1][0]
    )

    if i0 > 0:
        #qpair = (i0-1, abs(modes[i0-1][1]) + abs(modes[i0][1]))
        params["qWavePositon"] = zcfind(
            derivative,
            single=True,
            lb=modes[i0-1][0],
            rb=modes[i0][0]
        )

    i0 = maxpair[0]+1
    if i0+1 < len(modes):
        #spair = (i0+1, abs(modes[i0][1]) + abs(modes[i0+1][1]))
        params["sWavePositon"] = zcfind(
            derivative,
            single=True,
            lb=modes[i0][0],
            rb=modes[i0+1][0]
        )

    # Определение типа qrs-комплекса
    if params["qrsType"] is None:
        if params["rWavePositon"] is not None:
            if params["qWavePositon"] is not None:
                if params["sWavePositon"] is not None:
                    params["qrsType"] = "qRs"
                else:
                    params["qrsType"] = "qR"
            else:
                if params["sWavePositon"] is not None:
                    params["qrsType"] = "Rs"
                else:
                    params["qrsType"] = "R"

    return params, signcode


def find_points(x, fs, qrs_metadata, debug=True):

    bands = ddwt(x)
    modas = []
    for band in bands:
        moda = []
        # ищем положительные максимумы
        pos = argrelmax(band)[0]
        for i in pos:
            y = band[i]
            if y > 0:
                moda.append((i, y))
        # ищем отрицательные минимумы
        neg = argrelmin(band)[0]
        for i in neg:
            y = band[i]
            if y < 0:
                moda.append((i, y))

        moda.sort()
        modas.append(moda)

    # границы QRS здесь не определяем, надеемся на qrs_metadata

    summary = {}
    new_metadata = []
    noise = np.std(bands[1]) * 0.7
    print(noise)
    for qrs in qrs_metadata:

        pkdata = qrs.copy()

        r_scale = 2

        # начальный интервал для поиска волн
        lbound = int(qrs["qrs_start"] * fs)
        rbound = int(qrs["qrs_end"] * fs)

        modas_subset = filter(
            lambda x: lbound < x[0] < rbound and abs(x[1]) > noise,
            modas[r_scale]
        )

        params, codestr = pksearch(modas_subset, bands[r_scale])

        if debug:
            # для отладки: кодируем найденные фронты знаками +/-
            if codestr:
                summary[codestr] = summary.get(codestr, 0) + 1

        #pkdata["rwav"] = (lbound, rbound)
        pkdata.update(params)

        new_metadata.append(pkdata)

    print(json.dumps(summary, indent=1))
    return new_metadata