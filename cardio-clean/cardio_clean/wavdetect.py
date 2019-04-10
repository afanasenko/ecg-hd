# coding: utf-8

from scipy.signal import convolve, argrelmax, argrelmin, find_peaks
from metadata import *
from util import signal_channels
from matplotlib import pyplot as plt
from sigbind import detect_periodic


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

    return pk


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
            dly_lo = len(h)-1
            ap[:dly_lo] = ap[dly_lo] # хак
            # вместо прореживания сигналов (Маллат) расширяем характеристики
            # фильтров
            h = fexpand(h)
            g = fexpand(g)

    return approx, detail


def detect_all_extrema(modes, sig, smp_from, smp_to, marker):

    ex = []
    idx = -1

    for i, m in enumerate(modes[:-1]):

        x = m[0]
        y = m[1]
        xn = modes[i + 1][0]
        yn = modes[i + 1][1]

        if y*yn < 0:
            if y > yn:
                # find local max
                mpos = x + np.argmax(sig[x:xn])
                mval = sig[mpos]
                mpol = 1
            else:
                # find local min
                mpos = x + np.argmin(sig[x:xn])
                mval = sig[mpos]
                mpol = -1

            if mpos == marker:
                idx = len(ex)

            if smp_from <= mpos <= smp_to:
                # уберем двойные сразу
                if len(ex) and mpol == ex[-1][2]:
                    if mpol > 0:
                        if mval > ex[-1][1]:
                            ex[-1] = (mpos, mval, mpol)
                    else:
                        if mval < ex[-1][1]:
                            ex[-1] = (mpos, mval, mpol)
                else:
                    ex.append((mpos, mval, mpol))

    return ex, idx


def detect_r_pair(modes, smp_from, smp_to, bipol=False):

    maxd = 0      # максимальный перепад
    maxdpos = -1  # номер моды, дающей максимальный перепад
    i = 0
    for x, y in modes[:-1]:

        if smp_from <= x <= smp_to:
            if y * modes[i+1][1] < 0:
                diff = y - modes[i+1][1]

                if bipol:
                    diff = abs(diff)

                if diff > maxd:
                    maxd = diff
                    maxdpos = i
        i += 1

    return maxdpos, maxdpos+1


def qrssearch(modes, tight_bounds, approx, params, chan, isolevel,
              max_qrs_len):
    """

    :param modes: массив с модами производной
    :param tight_bounds: узкие границы для поиска центрального зубца
    :param approx: (сглаженный) входной сигнал
    :param detail:
    :param params: (inout) словарь с параметрами текущего qrs
    :param chan: номер канала
    :param isolevel: уровень изолинии для нахождения начала и конца R-зубца
    :param max_qrs_len: максимальное число отсчетов в QRS
    :return: none
    """

    # убираем все предыдущие значения для QRS
    # FIXME: хотя qrssearch не должен выполняться на уже размеченных сигналах
    erase_qrs(params, chan)

    if len(modes) < 2:
        return

    # Фронты самого крутого положительного зубца (ищем самую мощную пару +-)
    r0, r1 = detect_r_pair(modes, smp_from=tight_bounds[0],
                           smp_to=tight_bounds[1], bipol=False)

    # флаг показывает, есть ли у нас R-зубец
    have_r = r0 >= 0

    if not have_r:
        # не найдено ни одного положительного зубца, ищем отрицательный
        r0, r1 = detect_r_pair(modes, smp_from=tight_bounds[0],
                           smp_to=tight_bounds[1], bipol=True)

        if r0 < 0:
            return

        qs_from = modes[r0][0]
        qs_to = modes[r1][0]
        mpos = qs_from + np.argmin(approx[qs_from:qs_to])
        params["q_pos"][chan] = mpos
        params["r_pos"][chan] = None # тоже mpos ???
        params["s_pos"][chan] = mpos
        params["qrs_shape"][chan] = "qs"
        return

    # r0 - номер первой моды, соответствующей R-зубцу
    r_from = modes[r0][0]
    r_to = modes[r1][0]
    r_pos = r_from + np.argmax(approx[r_from:r_to])
    params["r_pos"][chan] = r_pos

    # уточненные границы комплекса
    qrs_from = r_pos - max_qrs_len/2
    qrs_to = r_pos + max_qrs_len/2

    e, r_idx = detect_all_extrema(modes, approx, qrs_from, qrs_to, r_pos)

    if r_idx < 0:  # ошибка в алгоритме, дальше ничего не выйдет все равно
        return

    q_pos = -1
    r1_pos = -1
    s1_pos = -1
    r2_pos = e[r_idx][0]
    s2_pos = -1

    num_peaks_right = len(e) - r_idx - 1

    # если 0 - значит нет правого S
    if num_peaks_right >= 1:
        if e[r_idx+1][2] < 0:
            s2_pos = e[r_idx+1][0]
        else:
            assert 0  # после перехода на find_peaks не должно быть двойных
            # однополярных пиков

    if num_peaks_right >= 3:
        r1_pos = e[r_idx][0]

        if e[r_idx+1][2] < 0:
            s1_pos = e[r_idx+1][0]
        else:
            assert 0  # после перехода на find_peaks не должно быть двойных
            # однополярных пиков

        if e[r_idx + 2][2] > 0:
            r2_pos = e[r_idx + 2][0]
        else:
            assert 0  # после перехода на find_peaks не должно быть двойных
            # однополярных пиков

        if e[r_idx + 3][2] < 0:
            s2_pos = e[r_idx + 3][0]
        else:
            assert 0  # после перехода на find_peaks не должно быть двойных
            # однополярных пиков

    num_peaks_left = r_idx
    # если 0 - значит нет Q
    if num_peaks_left >= 1:
        r1_pos = e[r_idx][0]
        if e[r_idx - 1][2] < 0:
            q_pos = e[r_idx - 1][0]
        else:
            assert 0  # после перехода на find_peaks не должно быть двойных
            # однополярных пиков

    if num_peaks_left >= 3 and num_peaks_right < 3:
        r2_pos = e[r_idx][0]

        if e[r_idx - 1][2] < 0:
            s1_pos = e[r_idx - 1][0]
        else:
            assert 0  # после перехода на find_peaks не должно быть двойных
            # однополярных пиков

        if e[r_idx - 2][2] > 0:
            r1_pos = e[r_idx - 2][0]
        else:
            assert 0

        if e[r_idx - 3][2] < 0:
            q_pos = e[r_idx - 3][0]
        else:
            assert 0

    if q_pos >= 0:
        params["q_pos"][chan] = q_pos

    if r1_pos == r2_pos:
        r2_pos = -1

    if r1_pos >= 0 and r2_pos >= 0:
        rm1 = approx[r1_pos] - isolevel
        rm2 = approx[r2_pos] - isolevel
        rm = max(rm1, rm2)
        t = 0.16
        if rm1 < t * rm:
            r1_pos = -1
        elif rm2 < t * rm:
            r2_pos = -1

    if r1_pos >= 0 and r2_pos >= 0:
        params["r_pos"][chan] = min(r1_pos, r2_pos)
        params["r2_pos"][chan] = max(r1_pos, r2_pos)
    elif r1_pos >= 0:
        params["r_pos"][chan] = r1_pos
    elif r2_pos >= 0:
        params["r_pos"][chan] = r2_pos
    else:
        assert 0

    if r2_pos < 0:
        if s1_pos < r_pos:
            s1_pos = -1
        elif s2_pos < r_pos:
            s2_pos = -1

    rmpos = max(r1_pos, r2_pos)
    if s1_pos > rmpos and s2_pos > rmpos:
        s1_pos = min(s1_pos, s2_pos)
        s2_pos = -1

    if s1_pos >= 0 and s2_pos >= 0:

        first_s = min(s1_pos, s2_pos)
        second_s = max(s1_pos, s2_pos)

        params["s_pos"][chan] = second_s
        params["s2_pos"][chan] = first_s

        # требует уточнения.- должен ли быть основной S всегда глубже
        # удаляем r2 и s2
        if approx[second_s] > approx[first_s]:
            params["s_pos"][chan] = first_s
            params["s2_pos"][chan] = None
            params["r2_pos"][chan] = None

    elif s1_pos >= 0:
        params["s_pos"][chan] = s1_pos
    elif s2_pos >= 0:
        params["s_pos"][chan] = s2_pos

    rpos = params["r_pos"][chan]
    if rpos is not None:
        r_thresh = 0.1 * abs(approx[rpos] - isolevel)
        x = rpos
        while abs(approx[x] - isolevel) > r_thresh:
            if x <= qrs_from:
                break
            x -= 1

        params["r_start"][chan] = x
        x = rpos
        while abs(approx[x] - isolevel) > r_thresh:
            if x >= qrs_to:
                break
            x += 1

        params["r_end"][chan] = x
    else:
        params["r_start"][chan] = None
        params["r_end"][chan] = None


def ptsearch(modes, approx, bias, limits, height):
    """
    Поиск зубцов P и T
    :param modes: список экстремумов производной
    :param approx: опорный сигнал для поиска пиков
    :param height: минимальное отклонение от изолинии
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

    wave_center, ampl = modefind(
        approx, lb=modes[i0][0], rb=modes[i0 + 1][0], bias=bias
    )

    if ampl < height:
        return None, None, None


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
    pos = start_idx + find_peaks(band[start_idx:end_idx+1], height=thresh)[0]
    for i in pos:
        moda.append((i, band[i]))
    # ищем отрицательные минимумы
    neg = start_idx + find_peaks(-band[start_idx:end_idx+1], height=thresh)[0]
    for i in neg:
        moda.append((i, band[i]))

    moda.sort()
    return moda


def pma_search(modes, smp_from, smp_to, max_dur):
    if len(modes) < 2:
        return 0

    # Фронты самого крутого положительного зубца (ищем самую мощную пару +-)
    r0, r1 = detect_r_pair(modes, smp_from=smp_from,
                           smp_to=smp_to, bipol=False)

    # флаг показывает, есть ли у нас R-зубец
    have_r = r0 >= 0

    if not have_r:
        # не найдено ни одного положительного зубца, ищем отрицательный
        r0, r1 = detect_r_pair(modes, smp_from=smp_from,
                           smp_to=smp_to, bipol=True)

        if r0 < 0:
            return 0

    dur = modes[r1][0] - modes[r0][0]
    if dur < max_dur:
        return int(0.5*(modes[r1][0] + modes[r0][0]))


def find_points(
        sig,
        fs,
        metadata,
        bias,
        gain,
        **kwargs):
    """
        Поиск характерных точек
    :param sig: входной сигнал (многоканальный)
    :param fs: частота дискретизации
    :param metadata: первичная сегментация, содержащая qrs_start, qrs_end,
    qrs_center
    :param bias: уровень изолинии
    :param gain: усиление
    :return: None (результатом являются измененные значения в metadata)
    """

    qrs_duration_max = int(fs*kwargs.get(
        "qrs_duration_max",
        config.WAVES["qrs_duration_max"]
    ))

    pma_detection_on = kwargs.get(
        "pma_detection_on",
        config.PACEMAKER["detection"]
    )

    pma_duration_max = fs*kwargs.get(
        "pma_duration_max",
        config.PACEMAKER["spike_duration_max"]
    )

    pma_scale = 1   # номер уровня для поиска артефактов кардиостимулятора
    r_scale = 2     # номер уровня для поиска R-зубца
    p_scale = 4     # номер уровня для поиска P-зубца
    t_scale = 4     # номер уровня для поиска T-зубца
    f_scale = 4     # номер уровня для обнаружения трепетаний

    t_window_fraction = 0.6
    p_window_fraction = 1.0 - t_window_fraction

    num_scales = max(r_scale, p_scale, t_scale)
    num_cycles = len(metadata)
    pilot_chan = 1 if sig.shape[1] > 1 else 0

    for chan, x in signal_channels(sig):

        approx, detail = ddwt(x-bias[chan], num_scales=num_scales)

        # границы QRS здесь не определяем, надеемся на metadata

        # очень приближенная оценка шума
        noise = np.std(detail[1]) * 0.7
        #print(noise)

        if chan == pilot_chan:
            fibpos = find_peaks(detail[f_scale], height=noise/2)[0]
        else:
            fibpos = []

        for ncycle, qrs in enumerate(metadata):

            prev_r = int(metadata[ncycle - 1]["qrs_center"] * fs) \
                if ncycle else 0
            next_r = int(metadata[ncycle + 1]["qrs_center"] * fs) \
                if ncycle < num_cycles - 1 else len(x)
            cur_r = int(qrs["qrs_center"] * fs)

            # оценка изолинии
            iso = np.percentile(approx[r_scale][prev_r:next_r], 15)
            qrs["isolevel"][chan] = iso / gain[chan]

            # Поиск зубцов Q, R, S
            # узкое окно для поиска только R
            this_qrs = [int(qrs["qrs_start"] * fs), int(qrs["qrs_end"] * fs)]

            tight_bounds = [
                max(0, this_qrs[0]),
                min(len(x) - 1, this_qrs[1])
            ]

            # более широкое окно для поиска остальных зубцов
            prev_qrs =\
                int(metadata[ncycle - 1]["qrs_end"] * fs) if ncycle else 0

            next_qrs =\
                int(metadata[ncycle + 1]["qrs_start"] * fs) \
                    if ncycle < num_cycles - 1 else len(x)

            loose_bounds = [
                int((tight_bounds[0] + prev_qrs)/2),
                int((tight_bounds[1] + next_qrs)/2)
            ]

            #if ncycle==34 and chan==0:
            #    fig, axarr = plt.subplots(2, 1, sharex="col")
            #    fleft = int(metadata[ncycle-1]["qrs_end"]*fs)
            #    fright = int(qrs["qrs_start"]*fs)
            #    #x1 = tight_bounds[0]
            #    #x2 = tight_bounds[1]
            #    x1 = fleft
            #    x2 = fright
            #    xval = np.arange(x1, x2)
            #    axarr[0].plot(xval, approx[r_scale][x1:x2])
            #    axarr[0].grid()
            #    acf, v = detect_periodic(detail[f_scale][x1:x2])
            #    xv = x1 + np.arange(0, len(acf))
            #    axarr[1].plot(xv, acf)
            #    #axarr[1].plot(xval, approx[f_scale][x1:x2])
            #    #axarr[1].plot(xval, detail[f_scale][x1:x2], "m")
            #    axarr[1].grid()
            #    print("Look at the plots")
            #    plt.show(block=False)
            #print(ncycle, chan)

            if pma_detection_on:
                pma_modes = find_extrema(
                    detail[pma_scale], loose_bounds[0], loose_bounds[1],
                    noise / 2
                )

                pma = pma_search(pma_modes, tight_bounds[0], tight_bounds[1],
                           max_dur=pma_duration_max)

                if pma:
                    qrs["pma"][chan].append(pma)

            # все пики производной в широком окне
            modes = find_extrema(
                detail[r_scale], loose_bounds[0], loose_bounds[1], noise/2
            )

            qrssearch(modes, tight_bounds, approx[r_scale],
                      qrs, chan, iso, qrs_duration_max)

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
                limits=pwindow,
                height=noise/2
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
                    qrs["isolevel"][chan] = iso / gain[chan]

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
                limits=twindow,
                height=noise
            )

            qrs["t_pos"][chan] = tcenter
            qrs["t_start"][chan] = tleft
            qrs["t_end"][chan] = tright

            # поиск F-волн в промежутках между qrs
            if chan == pilot_chan and ncycle:
                #fleft = int(metadata[ncycle-1]["qrs_end"]*fs)
                #fright = int(qrs["qrs_start"]*fs)

                fleft = get_cycle_end(metadata[ncycle-1], chan, fs)
                fright = get_cycle_start(qrs, chan, fs)


                fpeaks = [fpk for fpk in fibpos if
                                            fleft<fpk<fright]
                numf = len(fpeaks)
                qrs["f_waves"][chan] = numf

                # Берем промежуток между QRS (с возможным захватом P и T)
                # и обнаруживаем в нем периодичность
                rest_range = [
                    int(metadata[ncycle-1]["qrs_end"]*fs),
                    int(qrs["qrs_start"]*fs)
                    ]

                # защита от слишком коротких пауз
                # TODO: разобраться, почему это происходит
                if rest_range[1] - rest_range[0] > 8:
                    pk = detect_periodic(detail[f_scale][
                                         rest_range[0]:rest_range[1]])[1]
                else:
                    pk = 0

                qrs["flutter"][chan] = pk

