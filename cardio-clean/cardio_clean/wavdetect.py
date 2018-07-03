# coding: utf-8

import json
import numpy as np
from scipy.signal import convolve, argrelmax, argrelmin
from metadata import *


def modefind(x, lb=0, rb=0, bias=0.0):

    if rb:
        rb = min(rb, len(x))
    else:
        rb = len(x)

    pk = (0, 0.0)
    for i in range(lb, rb):
        val = abs(x[i]-bias)
        if val > pk[1]:
            pk = (i, val)

    return pk[0]


def zcfind(x, lb=0, rb=0):

    if rb:
        rb = min(rb, len(x))
    else:
        rb = len(x)

    for i in range(lb+1, rb):
        if x[i-1]*x[i] < 0:
            return i-1


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
    Растягивание массива вдвое с добавлением нулей между элементами (
    используется в ДВП)
    :param x: одномерный массив из N элементов
    :return: расширенный массив из 2N элементов
    """
    n = len(x)
    y = np.array([0]*2*n, float)
    for i in range(len(y)):
        if i%2 == 0:
            y[i] = x[i/2]

    return y


def ddwt(x, num_scales):
    """
    Дискретное вейвлет-преобразование без прореживания
    :param x: входной сигнал
    :param num_scales: число уровней разложения
    :return:
    """
    h = np.array([1, 3, 3, 1], float) / 8
    g = np.array([2, -2], float)
    signal_len = len(x)

    detail = []
    approx = []
    ap = x.copy()
    detail.append(ap.copy())  # на нулевом уровне храним исходный сигнал
    approx.append([])

    for s in range(num_scales):
        dly = 2**s
        hif = convolve(ap, g, mode="full")[dly:dly+signal_len]
        detail.append(hif)
        approx.append(ap)
        if s < num_scales-1:
            ap = convolve(ap, h, mode="full")[dly:dly+signal_len]
            # вместо прореживания сигналов (Маллат) расширяем характеристики
            # фильтров
            h = fexpand(h)
            g = fexpand(g)

    return approx, detail


def qrssearch(modes, derivative, params):
    """

    :param modes:
    :param derivative:
    :param params: (inout) словарь с параметрами текущего qrs
    :return: символический код
    """

    signcode = ""

    if len(modes) < 2:
        return signcode

    # кодируем найденные фронты знаками +/-
    signcode = "".join(["+" if x[1] > 0 else "-" for x in modes])

    if signcode == "-+":
        params["qrsType"] = "qs"

    # Фронты самого мощного (R или  qs) зубца
    maxpair = (0,0)
    for i, posval in enumerate(modes):
        if i and posval[1]*modes[i - 1][1] < 0:
            diff = abs(posval[1]) + abs(modes[i - 1][1])
            if diff > maxpair[1]:
                maxpair = (i-1, diff)

    i0 = maxpair[0]

    params["waves"]["r"]["center"] = zcfind(
        derivative,
        lb=modes[i0][0],
        rb=modes[i0 + 1][0]
    )

    if i0 > 0:
        params["waves"]["q"]["center"] = zcfind(
            derivative,
            lb=modes[i0-1][0],
            rb=modes[i0][0]
        )

    i0 = maxpair[0]+1
    if i0+1 < len(modes):
        params["waves"]["s"]["center"] = zcfind(
            derivative,
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

    return signcode


def edgefind(x0, y0, dx, dy, bias):
    return int(x0 - dx * (y0-bias) / dy)


def ptsearch(modes, approx, bias=0.0):
    """
    Поиск зубцов P и T
    :param modes: список экстремумов производной
    :param approx: опорный сигнал для поиска пиков
    :return: wave_left, wave_center, wave_right
    """

    if len(modes) < 2:
        return None, None, None

    # Фронты самого мощного зубца
    maxpair = (0,0)
    for i, posval in enumerate(modes):
        if i and posval[1]*modes[i - 1][1] < 0:
            diff = abs(posval[1]) + abs(modes[i - 1][1])
            if diff > maxpair[1]:
                maxpair = (i-1, diff)

    i0 = maxpair[0]

    wave_center = modefind(
        approx, lb=modes[i0][0], rb=modes[i0 + 1][0], bias=bias
    )

    # строим касательную в наиболее крутой точке переднего фронта

    x0 = modes[i0][0]
    y0 = approx[x0]
    dy = approx[x0+1] - approx[x0-1]

    if dy > 0:
        wave_left = edgefind(x0, y0, 2.0, dy, bias)
        if wave_left >= wave_center or wave_left <= 0:
            wave_left = None
    else:
        wave_left = None

    # строим касательную в наиболее крутой точке заднего фронта

    x0 = modes[i0+1][0]
    y0 = approx[x0]
    dy = approx[x0+1] - approx[x0-1]

    if dy < 0:
        wave_right = edgefind(x0, y0, 2.0, dy, bias)
        if wave_right <= wave_center or wave_right >= len(approx):
            wave_right = None
    else:
        wave_right = None

    return wave_left, wave_center, wave_right


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


def find_points(
        x,
        fs,
        qrs_metadata,
        bias=0.0,
        debug=False):
    """
        Поиск характерных точек
    :param x: входной сигнал (1 канал)
    :param fs: частота дискретизации
    :param qrs_metadata: первичная сегментация: qrs_start, qrs_end, r_position
    :param bias: уровень изолинии
    :param debug:
    :return: список словарей метаданных по каждому циклу
    """

    r_scale = 2
    p_scale = 4
    t_scale = 4
    t_window_fraction = 0.6
    p_window_fraction = 1.0 - t_window_fraction

    approx, detail = ddwt(x-bias, num_scales=max(r_scale, p_scale,
                                                     t_scale))
    modas = []
    for band in detail:
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
    metadata = []

    # очень приближенная оценка шума
    noise = np.std(detail[1]) * 0.7
    if debug:
        print(noise)

    for ncycle, qrs in enumerate(qrs_metadata):

        pkdata = metadata_new()

        # Поиск зубцов Q, R, S
        # окно для поиска
        lbound = int(qrs["qrs_start"] * fs)
        rbound = int(qrs["qrs_end"] * fs)

        pkdata["qrs_start"] = qrs["qrs_start"]
        pkdata["qrs_end"] = qrs["qrs_end"]

        modas_subset = range_filter( modas[r_scale], lbound, rbound, noise/2)

        codestr = qrssearch(modas_subset, detail[r_scale], pkdata)

        if debug:
            # для отладки: кодируем найденные фронты знаками +/-
            if codestr:
                summary[codestr] = summary.get(codestr, 0) + 1

        prev_r = int(qrs_metadata[ncycle - 1]["r_position"][0] * fs) \
            if ncycle else 0
        next_r = int(qrs_metadata[ncycle + 1]["r_position"][0] * fs) \
            if ncycle < len(qrs_metadata) - 1 else len(x)
        cur_r = int(qrs["r_position"][0] * fs)

        # оценка изолинии
        iso = np.percentile(approx[r_scale][prev_r:next_r], 15)
        pkdata["isolevel"] = iso

        # поиск P-зубца
        # окно для поиска
        wlen = (cur_r - prev_r) * p_window_fraction
        pwindow = [
            int(prev_r + wlen),
            cur_r
        ]

        modas_subset = range_filter(
            modas[p_scale], pwindow[0], pwindow[1], noise/2
        )

        # последняя мода перед R не учитывается, потому что относится к QRS
        pleft, pcenter, pright = ptsearch(
            modas_subset[:-1],
            approx[r_scale+1],
            bias=iso
        )

        pkdata["waves"]["p"]["center"] = pcenter
        pkdata["waves"]["p"]["start"] = pleft
        pkdata["waves"]["p"]["end"] = pright

        # поиск T-зубца
        # окно для поиска
        wlen = (next_r - cur_r) * t_window_fraction
        twindow = [
            cur_r,
            int(cur_r + wlen)
        ]

        modas_subset = range_filter(
            modas[p_scale], twindow[0], twindow[1], noise / 4
        )

        # первая мода справа от R не учитывается, потому что относится к QRS
        tleft, tcenter, tright = ptsearch(
            modas_subset[1:],
            approx[r_scale+1],
            bias=iso
        )

        pkdata["waves"]["t"]["center"] = tcenter
        pkdata["waves"]["t"]["start"] = tleft
        pkdata["waves"]["t"]["end"] = tright

        #print(json.dumps(pkdata, indent=1))
        metadata.append(pkdata)

    if debug:
        print(json.dumps(summary, indent=1))

    return metadata