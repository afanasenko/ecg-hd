# coding: utf-8

from scipy.signal import convolve, argrelmax, argrelmin
from metadata import *
from util import signal_channels
from matplotlib import pyplot as plt

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
            dly_lo = len(h)-1
            ap = convolve(ap, h, mode="full")[dly_lo:dly_lo+signal_len]
            # вместо прореживания сигналов (Маллат) расширяем характеристики
            # фильтров
            h = fexpand(h)
            g = fexpand(g)

    return approx, detail


def qrssearch(modes, approx, derivative, params, chan, isolevel):
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
    lb = modes[0][0]
    rb = modes[-1][0]
    for i, posval in enumerate(modes):
        if i and posval[1]*modes[i - 1][1] < 0:
            diff = abs(posval[1]) + abs(modes[i - 1][1])
            if diff > maxpair[1]:
                maxpair = (i-1, diff)

    i0 = maxpair[0]

    rpos = zcfind(
        derivative,
        lb=modes[i0][0],
        rb=modes[i0 + 1][0]
    )
    params["r_pos"][chan] = rpos

    if rpos is not None:
        r_thresh = 0.1 * abs(approx[rpos] - isolevel)
        x = rpos
        while abs(approx[x] - isolevel) > r_thresh:
            if x <= lb:
                break
            x -= 1

        params["r_start"][chan] = x
        x = rpos
        while abs(approx[x] - isolevel) > r_thresh:
            if x >= rb:
                break
            x += 1

        params["r_end"][chan] = x
    else:
        params["r_start"][chan] = None
        params["r_end"][chan] = None


    q_search_rb = modes[i0][0]
    if params["r_start"][chan] is not None:
        q_search_rb = min(q_search_rb, params["r_start"][chan])

    if i0 > 0:
        params["q_pos"][chan] = zcfind(
            derivative,
            lb=modes[i0-1][0],
            rb=q_search_rb
        )

    i0 = maxpair[0]+1

    s_search_lb = modes[i0][0]
    if params["r_end"][chan] is not None:
        s_search_lb = min(s_search_lb, params["r_end"][chan])

    if i0+1 < len(modes):
        params["s_pos"][chan] = zcfind(
            derivative,
            lb=s_search_lb,
            rb=modes[i0+1][0]
        )

    # Определение типа qrs-комплекса по II-му стандартному отведению
    if chan == 1:
        if params["qrsType"] is None:
            if params["r_pos"][chan] is not None:
                if params["q_pos"][chan] is not None:
                    if params["s_pos"][chan] is not None:
                        params["qrsType"] = "qRs"
                    else:
                        params["qrsType"] = "qR"
                else:
                    if params["s_pos"][chan] is not None:
                        params["qrsType"] = "Rs"
                    else:
                        params["qrsType"] = "R"

    return signcode


def ptsearch(modes, approx, bias, limits):
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
    y0 = approx[x0] - bias
    dy = approx[x0+1] - approx[x0-1]

    if abs(approx[x0+1] - bias) > abs(approx[x0-1] - bias):
        wave_left = int(x0 - 2.0 * y0 / dy)
        if wave_left >= wave_center or wave_left <= limits[0]:
            wave_left = None
    else:
        wave_left = None

    # строим касательную в наиболее крутой точке заднего фронта

    x0 = modes[i0+1][0]
    y0 = approx[x0] - bias
    dy = approx[x0+1] - approx[x0-1]

    if abs(approx[x0+1] - bias) < abs(approx[x0-1] - bias):
        wave_right = int(x0 - 2.0 * y0 / dy)
        if wave_right <= wave_center or wave_right >= limits[1]:
            wave_right = None
    else:
        wave_right = None

    return wave_left, wave_center, wave_right


def find_extrema(band, start_idx, end_idx, thresh):

    moda = []
    # ищем положительные максимумы
    pos = start_idx + argrelmax(band[start_idx:end_idx+1])[0]
    for i in pos:
        y = band[i]
        if y > thresh:
            moda.append((i, y))
    # ищем отрицательные минимумы
    neg = start_idx + argrelmin(band[start_idx:end_idx+1])[0]
    for i in neg:
        y = band[i]
        if y < -thresh:
            moda.append((i, y))

    moda.sort()
    return moda


def find_points(
        sig,
        fs,
        metadata,
        bias,
        debug=False):
    """
        Поиск характерных точек
    :param sig: входной сигнал (многоканальный)
    :param fs: частота дискретизации
    :param metadata: первичная сегментация, содержащая qrs_start, qrs_end,
    qrs_center
    :param bias: уровень изолинии
    :param debug: отладочный флаг
    :return: None (результатом являются измененные значения в metadata)
    """

    r_scale = 2
    p_scale = 4
    t_scale = 4
    t_window_fraction = 0.6
    p_window_fraction = 1.0 - t_window_fraction

    num_scales = max(r_scale, p_scale, t_scale)

    for chan, x in signal_channels(sig):

        approx, detail = ddwt(x-bias[chan], num_scales=num_scales)

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

            moda.sort(key=lambda x: x[0])
            modas.append(moda)

        # для обнаружения трепетаний
        fscale = 4
        if chan == 1:
            fibpos = argrelmin(detail[fscale], order=3)[0]
        else:
            fibpos = []

        # границы QRS здесь не определяем, надеемся на metadata

        #plt.plot(x)
        #plt.show(block=False)

        # очень приближенная оценка шума
        noise = np.std(detail[1]) * 0.7
        if debug:
            print(noise)

        for ncycle, qrs in enumerate(metadata):

            prev_r = int(metadata[ncycle - 1]["qrs_center"] * fs) \
                if ncycle else 0
            next_r = int(metadata[ncycle + 1]["qrs_center"] * fs) \
                if ncycle < len(metadata) - 1 else len(x)
            cur_r = int(qrs["qrs_center"] * fs)

            # оценка изолинии
            #iso = np.median(approx[r_scale][prev_r:next_r])
            iso = np.percentile(approx[r_scale][prev_r:next_r], 15)
            qrs["isolevel"][chan] = iso

            # Поиск зубцов Q, R, S
            # окно для поиска
            addsmp = 3
            lbound = max(0, int(qrs["qrs_start"] * fs) - addsmp)
            rbound = min(sig.shape[0]-1, int(qrs["qrs_end"] * fs) + addsmp)

            modas_subset = find_extrema(
                detail[r_scale], lbound, rbound, noise/2
            )

            qrssearch(modas_subset, approx[r_scale], detail[r_scale], qrs,
                      chan, iso)

            # поиск P-зубца
            # окно для поиска
            wlen = (cur_r - prev_r) * p_window_fraction

            p_search_lb = int(prev_r + wlen)
            if ncycle:
                prev_t = metadata[ncycle-1]["t_end"][chan]
                if prev_t is None:
                    prev_t = metadata[ncycle - 1]["t_pos"][chan]
                if prev_t is not None:
                    p_search_lb = max(p_search_lb,  prev_t)

            pwindow = [
                p_search_lb,
                cur_r
            ]

            modas_subset = find_extrema(
                detail[p_scale], pwindow[0], pwindow[1], noise/2
            )

            # последняя мода перед R не учитывается, потому что относится к QRS
            pleft, pcenter, pright = ptsearch(
                modas_subset[:-1],
                approx[r_scale+1],
                bias=iso,
                limits=pwindow
            )

            qrs["p_pos"][chan] = pcenter
            qrs["p_start"][chan] = pleft
            qrs["p_end"][chan] = pright

            # уточнение уровня изолинии по интервалу PQ
            if pright is not None:
                pq_end = qrs["q_pos"][chan]
                if pq_end is None:
                    pq_end = qrs["r_start"][chan]
                if pq_end is None:
                    pq_end = qrs["r_pos"][chan]
                if pq_end is not None and pq_end - pright > 1:
                    iso = np.median(approx[r_scale][pright:pq_end])
                    qrs["isolevel"][chan] = iso

            # поиск T-зубца
            # окно для поиска
            wlen = (next_r - cur_r) * t_window_fraction
            twindow = [
                cur_r,
                int(cur_r + wlen)
            ]

            modas_subset = find_extrema(
                detail[p_scale], twindow[0], twindow[1], noise / 4
            )

            # первая мода справа от R не учитывается, потому что относится к QRS
            tleft, tcenter, tright = ptsearch(
                modas_subset[1:],
                approx[r_scale+1],
                bias=iso,
                limits=twindow
            )

            qrs["t_pos"][chan] = tcenter
            qrs["t_start"][chan] = tleft
            qrs["t_end"][chan] = tright

            # поиск F-волн в промежутках между qrs
            if ncycle:
                fleft = int(metadata[ncycle-1]["qrs_end"]*fs)
                fright = int(qrs["qrs_start"]*fs)
                numf = 0

                for f in fibpos:
                    if fleft <= f <= fright:
                        if detail[fscale][f] < - noise:
                            numf += 1
                    elif f > fright:
                        break
                qrs["f_waves"][chan] = numf

                #if numf > 1:
                #   plt.plot(approx[1]-iso)
                #   plt.plot(detail[fscale] + 1)
                #   plt.xlim((fleft-fs,fright+fs))
                #   plt.plot([fleft, fright], [-2*noise, -2*noise])
                #   plt.show()


def detect_f_waves(x):

    pass

