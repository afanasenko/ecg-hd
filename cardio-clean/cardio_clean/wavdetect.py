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
        return zc[0] if zc else None
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
    Растягивание массива вдвое с добавлением нулей между элементами
    :param x: одномерный массив из N элементов
    :return: расширенный массив из 2N элементов
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
    signal_len = len(x)

    decomposition = []
    ap = x.copy()
    decomposition.append(ap.copy())  # на нулевом уровне храним исходный сигнал

    for s in range(num_scales):
        dly = 2**s
        hif = convolve(ap, g, mode="full")[dly:dly+signal_len]
        decomposition.append(hif)
        if s < num_scales-1:
            ap = convolve(ap, h, mode="full")[dly:dly+signal_len]
            # вместо прореживания сигналов (Маллат) расширяем характеристики
            # фильтров
            h = fexpand(h)
            g = fexpand(g)

    return decomposition


def pksearch(modes, derivative):

    params = {
        "waves": {},
        "qrsType": None,
    }
    params["waves"].update(makewave("q"))
    params["waves"].update(makewave("r"))
    params["waves"].update(makewave("s"))

    signcode = ""

    if len(modes) < 2:
        return params, signcode

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

    params["waves"]["r"]["center"] = zcfind(
        derivative,
        single=True,
        lb=modes[i0][0],
        rb=modes[i0 + 1][0]
    )

    if i0 > 0:
        params["waves"]["q"]["center"] = zcfind(
            derivative,
            single=True,
            lb=modes[i0-1][0],
            rb=modes[i0][0]
        )

    i0 = maxpair[0]+1
    if i0+1 < len(modes):
        params["waves"]["s"]["center"] = zcfind(
            derivative,
            single=True,
            lb=modes[i0][0],
            rb=modes[i0+1][0]
        )

    # Определение типа qrs-комплекса
    if params["qrsType"] is None:
        if params["waves"]["r"]["center"] is not None:
            if params["waves"]["q"]["center"] is not None:
                if params["waves"]["s"]["center"] is not None:
                    params["qrsType"] = "qRs"
                else:
                    params["qrsType"] = "qR"
            else:
                if params["waves"]["s"]["center"] is not None:
                    params["qrsType"] = "Rs"
                else:
                    params["qrsType"] = "R"

    return params, signcode


def makewave(name, pos=None, height=None, start=None, end=None):
    return {name: {
        "start": start, "end": end, "center": pos, "height": height
    }}


def ptsearch(modes, derivative):
    """
    Поиск зубцов P и T
    :param modes:
    :param derivative:
    :return:
    """

    if len(modes) < 2:
        return None

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


    # Фронты самого мощного зубца
    maxpair = (0,0)
    for i, posval in enumerate(modes):
        if i:
            diff = abs(posval[1]) + abs(modes[i - 1][1])
            if diff > maxpair[1]:
                maxpair = (i-1, diff)

    i0 = maxpair[0]

    p_wave_center = zcfind(
        derivative,
        single=True,
        lb=modes[i0][0],
        rb=modes[i0 + 1][0]
    )

    return p_wave_center


def range_filter(x, lb, rb, thresh):
    ret = []
    for m0, m1 in x:
        if m0 < lb:
            continue
        if m0 > rb:
            break

        if abs(m1) > thresh:
            ret.append((m0, m1))

    return ret


def find_points(x, fs, qrs_metadata, j_offset_ms=60, debug=False):

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

    # очень приближенная оценка шума
    noise = np.std(bands[1]) * 0.7
    if debug:
        print(noise)

    for ncycle, qrs in enumerate(qrs_metadata):

        pkdata = qrs.copy()

        # Поиск зубцов Q, R, S
        r_scale = 2
        # окно для поиска
        lbound = int(qrs["qrs_start"] * fs)
        rbound = int(qrs["qrs_end"] * fs)

        modas_subset = range_filter( modas[r_scale], lbound, rbound, noise)

        params, codestr = pksearch(modas_subset, bands[r_scale])
        pkdata.update(params)

        if debug:
            # для отладки: кодируем найденные фронты знаками +/-
            if codestr:
                summary[codestr] = summary.get(codestr, 0) + 1

        # поиск P-зубца

        p_scale = 4
        # окно для поиска
        prev_r = int(qrs_metadata[ncycle-1]["r_wave_center"][0] * fs) \
            if ncycle else 0
        cur_r = int(qrs["r_wave_center"][0] * fs)
        pwindow = [
            int((prev_r + cur_r)/2),
            cur_r
        ]

        modas_subset = range_filter(modas[p_scale], pwindow[0], pwindow[1],
                                    noise/2)

        # последняя мода не учитывается, потому что относится к QRS
        p_wave_center = ptsearch(modas_subset[:-1], bands[p_scale])

        pkdata["waves"].update(makewave("p", p_wave_center))

        # поиск T-зубца

        t_scale = 4
        # окно для поиска
        next_r = int(qrs_metadata[ncycle+1]["r_wave_center"][0] * fs) \
            if ncycle < len(qrs_metadata)-1 else len(x)

        cur_r = int(qrs["r_wave_center"][0] * fs)
        twindow = [
            cur_r,
            int((next_r + cur_r)/2)
        ]

        modas_subset = range_filter(modas[p_scale], twindow[0], twindow[1],
                                    noise / 2)

        # первая мода не учитывается, потому что относится к QRS
        t_wave_center = ptsearch(modas_subset[1:], bands[t_scale])

        pkdata["waves"].update(makewave("t", t_wave_center))

        # точка J не обнаруживается, а ставится со смещением от R-зубца
        rc = pkdata["waves"]["r"]["center"]
        if rc is not None:
            j_point = rc + int(fs*j_offset_ms/1000.0)
            if j_point > len(x) - 1:
                j_point = None
        else:
            j_point = None

        pkdata["waves"].update(makewave("j", j_point))

        # запись высоты зубцов
        for wave in pkdata["waves"]:
            pos = pkdata["waves"][wave]["center"]
            if pos is not None:
                pkdata["waves"][wave]["height"] = x[pos]

        new_metadata.append(pkdata)

    if debug:
        print(json.dumps(summary, indent=1))

    return new_metadata