# coding: utf-8

from util import signal_channels
from scipy.stats import linregress

"""
Метаданные подразделяются на первичные и вторичные.
Первичные являются результатом автоматической сегментации и могут быть
скорректированы вручную. Для расчета вторичных данных используются как
первичные, так и сам сигнал. После ручного редактирования необходимо
пересчитывать вторичные данные, не проводя повторной сегментации сигнала.

# ######################################
Первичные метаданные включают:

Q, R, S : центр зубца (метки времени)
P, Т: начало, центр, конец (метки времени)

qrsType: вид комплекса. Обнаруживаемые разновидности qrs - R, qR, Rs, qs, qRs.


Артефакты кодируются значением qrsType==None

# ######################################

Вторичные данные по зубцам
P, Q, R, S, T (амплитуда)
T: крутизна переднего фронта, крутизна заднего фронта, симметрия, острота

Вторичные данные по ритму
RR-интервал в миллисекундах
ЧСС (число ударов в минуту)

Вторичные данные по ST-сегменту
начало ST-сегмента (точка J): положение, уровень
точка ST+ (J + 0,08с): положение, уровень
наклон сегмента
смещение сегмента от изолинии
длительность сегмента ST в миллисекундах


Все временнЫе метки записываются как номер отсчета сигнала от начала
записи (тип int32)
Длительности и интервалы записываются в миллисекундах (тип float)
Все амплитудные параметры записываются в той же размерности, что и входной
сигнал (тип float).
"""

import numpy as np


def metadata_new(num_channels):
    return {

        # qrs-комплекс в целом
        "qrs_start": None,  # [секунд от начала записи] float
        "qrs_end": None,  # [секунд от начала записи] float
        "qrs_center": None,  # [секунд от начала записи] float
        "qrs_class_id": None,  # код класса string
        "qrsType": None,  # string
        "flags": "",  # string флаги ''(обычный)|'A'(артефакт)|'E'(
        # экстрасистола)|
        "complex_type": "U",

        # отдельные зубцы
        "p_start": [None]*num_channels,  # int array
        "p_end": [None]*num_channels,  # int array
        "p_pos": [None]*num_channels,  # int array
        "p_height": [None]*num_channels,  # float array
        "q_pos": [None]*num_channels,  # int array
        "q_height": [None]*num_channels,  # float array
        "r_start": [None]*num_channels,  # int array
        "r_end": [None]*num_channels,  # int array
        "r_pos": [None]*num_channels,  # int array
        "r_height": [None]*num_channels,  # float array
        "s_pos": [None]*num_channels,  # int array
        "s_height": [None]*num_channels,  # float array
        "t_start": [None]*num_channels,  # int array
        "t_end": [None]*num_channels,  # int array
        "t_pos": [None]*num_channels,  # int array
        "t_height": [None]*num_channels,  # float array

        # параметры ритма
        "RR": None,  # [миллисекунды] float
        "heartrate": None,  # [удары в минуту] float

        # оценка уровня изолинии
        "isolevel": [None]*num_channels,  # float array

        # ST-сегмент
        "st_start": [None]*num_channels,  # int array
        "st_plus": [None]*num_channels,  # int array
        "st_end": [None]*num_channels,  # int array
        "st_start_level": [None]*num_channels,  # float array
        "st_plus_level": [None]*num_channels,  # float array
        "st_end_level": [None]*num_channels,  # float array
        "st_offset": [None]*num_channels,  # float array
        "st_duration": [None]*num_channels,  # float array
        "st_slope": [None]*num_channels,  # float array

        # QT-интервал
        "qt_duration": [None]*num_channels,  # float array
        "qtc_duration": [None]*num_channels  # float array
    }


def samples_to_sec(smp, fs):
    return float(smp) / fs


def samples_to_ms(smp, fs):
    return smp * 1000.0 / fs


def ms_to_samples(ms, fs):
    return int(ms * fs / 1000.0)


def level_from_pos(d, chan, pos_key, val_key, sig, bias, gain):
    pos = d[pos_key][chan]
    if pos is None:
        d[val_key][chan] = None
    else:
        d[val_key][chan] = (sig[pos] - bias) / gain


def is_pvc(cycledata):
    """
    Проверка признака экстрасистолы в данном комплексе
    :param cycledata:
    :return: bool
    """
    return "E" in cycledata["flags"]


def is_artifact(cycledata):
    """
    Проверка признака артефакта в данном комплексе
    :param cycledata:
    :return: bool
    """
    return "A" in cycledata["flags"]


def set_artifact(cycledata):
    """
    Установка признака артефакта в данном комплексе
    :param cycledata:
    :return: bool
    """
    if "A" not in cycledata["flags"]:
        cycledata["flags"] += "A"


def reset_artifact(cycledata):
    """
    Сброс признака артефакта в данном комплексе
    :param cycledata:
    :return: bool
    """
    cycledata["flags"].replace("A", "")


def safe_r_pos(cycledata):
    heartbeat_channel=1

    # RR и ЧСС
    rz = cycledata["r_pos"][heartbeat_channel]
    if rz is None:
        realr = [x for x in cycledata["r_pos"] if x is not None]
        if len(realr) > 1:
            rz = np.median(realr)

    if rz is None:
        qz = cycledata["q_pos"][heartbeat_channel]
        sz = cycledata["s_pos"][heartbeat_channel]
        if qz is not None and sz is not None:
            rz = (qz + sz)/2
    return rz


def find_rr(metadata, pos, fs):

    cycledata = metadata[pos]
    num_cycles = len(metadata)

    # Не считаем RR для артефактов
    if is_artifact(cycledata):
        return

    # RR и ЧСС
    rz = safe_r_pos(cycledata)

    if rz is None:
        return

    vice = None

    if pos < num_cycles - 1:
        cand = pos + 1

        while cand < num_cycles:
            if is_artifact(metadata[cand]):
                continue
            vice = safe_r_pos(metadata[cand])
            if vice is not None:
                break
            cand += 1
    else:
        cand = pos - 1

        while cand > 0:
            if is_artifact(metadata[cand]):
                continue
            vice = safe_r_pos(metadata[cand])
            if vice is not None:
                break
            cand -= 1

    if vice is not None and vice != rz:
        return float(abs(rz - vice)) / fs


def metadata_postprocessing(metadata, sig, header, **kwargs):
    """
    Расчет вторичных параметров сигнала во всех отведениях

    Поскольку источник входных метаданных неизвестен, необходимо
    перезаписать значения всех зависимых ключей.
    :param metadata:
    :param sig:
    :param header: структура с полями fs, adc_gain, baseline
    :param kwargs: константы j_offset, jplus_offset_ms, min_st_ms
    :return: None (результатом являются измененные значения в metadata)
    """

    fs = header["fs"]

    j_offset_ms = kwargs.get("j_offset", 60)
    jplus_offset_ms = kwargs.get("jplus_offset", 80)

    numch = sig.shape[1] if sig.ndim == 2 else 1

    # ритм оценивается всегда по второму отведению
    heartbeat_channel = 1 if numch > 1 else 0

    # классификация тоже по второму отведению
    classification_channel = 1 if numch > 1 else 0

    num_cycles = len(metadata)

    # Удаление выбросов
    for ncycle, cycledata in enumerate(metadata):
        delta = ms_to_samples(60, fs)
        remove_outliers(cycledata, "p_pos", ("p_start", "p_end"), delta)
        remove_outliers(cycledata, "q_pos", [], delta)
        remove_outliers(cycledata, "r_pos", ("r_start", "r_end"), delta)
        remove_outliers(cycledata, "s_pos", [], delta)
        remove_outliers(cycledata, "t_pos", ("t_start", "t_end"), delta)

    for ncycle, cycledata in enumerate(metadata):

        # ######################################
        # RR и ЧСС
        rr = find_rr(metadata, ncycle, fs)
        if rr is None:
            set_artifact(cycledata)
            cycledata["RR"] = None
            cycledata["heartrate"] = None
        else:
            cycledata["RR"] = rr
            cycledata["heartrate"] = 60.0 / rr

        for chan, x in signal_channels(sig):

            # ######################################
            # точки J и J+
            # ставится со смещением от R-зубца
            rc = cycledata["r_pos"][chan]
            if rc is not None:
                j_point = rc + ms_to_samples(j_offset_ms, fs)
                # J не может быть раньше конца R
                r_end = cycledata["r_end"][chan]
                if r_end is not None:
                    j_point = max(j_point, r_end)
                if j_point > len(x) - 1:
                    j_point = None
                    jplus_point = None
                else:
                    jplus_point = j_point + ms_to_samples(jplus_offset_ms, fs)
                    if jplus_point > len(x) - 1:
                        jplus_point = None

                    elif cycledata["t_start"][chan] is not None:
                        tstart = cycledata["t_start"][chan]
                        jplus_point = min(j_point, tstart)

            else:
                j_point = None
                jplus_point = None

            # ######################################
            # ST
            st_end = cycledata["t_start"][chan]
            cycledata["st_start"][chan] = j_point
            cycledata["st_plus"][chan] = jplus_point
            cycledata["st_end"][chan] = st_end

            # ######################################
            # запись высоты зубцов
            bias = cycledata["isolevel"][chan]
            gain = header["adc_gain"][chan]

            level_from_pos(cycledata, chan, "p_pos", "p_height", x, bias, gain)
            level_from_pos(cycledata, chan, "q_pos", "q_height", x, bias, gain)
            level_from_pos(cycledata, chan, "r_pos", "r_height", x, bias, gain)
            level_from_pos(cycledata, chan, "s_pos", "s_height", x, bias, gain)
            level_from_pos(cycledata, chan, "t_pos", "t_height", x, bias, gain)
            level_from_pos(cycledata, chan, "st_start",
                           "st_start_level", x, bias, gain)
            level_from_pos(cycledata, chan, "st_plus", "st_plus_level", x,
                           bias, gain)
            level_from_pos(cycledata, chan, "st_end", "st_end_level", x,
                           bias, gain)

            # ######################################
            # ST (продолжение)
            if all((j_point, st_end)):
                dur = samples_to_ms(st_end - j_point, fs)
                if dur > kwargs.get("min_st_ms", 40):
                    cycledata["st_duration"][chan] = dur

                    cycledata["st_offset"][chan] = (np.mean(
                        sig[j_point:st_end]) - bias) / gain

                    cycledata["st_slope"][chan] = \
                        (cycledata["st_end_level"][chan] -
                         cycledata["st_start_level"][chan]) / \
                        cycledata["st_duration"][chan]
                else:
                    cycledata["st_duration"][chan] = None

            # ######################################
            # QT
            qt_start = cycledata["q_pos"][chan]
            if qt_start is None:
                qt_start = cycledata["r_start"][chan]

            qt_end = cycledata["t_end"][chan]

            if qt_start is not None and qt_end is not None:
                cycledata["qt_duration"][chan] = samples_to_sec(
                    qt_end - qt_start, fs)
            else:
                cycledata["qt_duration"][chan] = None

            # QTc

            qt = cycledata["qt_duration"][chan]
            if qt is not None and cycledata["RR"] is not None:
                cycledata["qtc_duration"][chan] = qt / np.sqrt(
                    cycledata["RR"])
            else:
                cycledata["qtc_duration"][chan] = None

            # ######################################
            #
            if chan == classification_channel:
                cycledata["complex_type"] = define_complex(
                    cycledata, fs, classification_channel
                )

    detect_pvc(metadata)


def is_sinus_cycle(meta):
    """
    Проверка на синусовый цикл
    :param meta:
    :return: True / False
    """

    std2 = 1

    # R-зубцы присутствуют во всех отведениях
    if all(meta["r_pos"]):
        for i, p in enumerate(meta["p_pos"]):
            # P-зубцы присутствуют во всех отведениях
            if p is None:
                return False
            # P-зубцы идут перед R-зубцами
            if meta["r_pos"][i] < p:
                return False

        # P-зубец во II-м отведении  положительный
        if meta["p_height"][std2] > 0:
            return True

        return False


def define_complex(meta, fs, channel):
    """

    :param meta:
    :param channel:
    :return: N - син. V - желуд. S - наджелуд. U - неизвестный
    """

    if is_sinus_cycle(meta):
        return "N"

    qrslen = estimate_qrslen(meta, fs, channel)

    # наджелудочковые комплексы - обычные, с P-зубцом
    # желудочковые комплексы - широкие, похожие на период синусоиды
    ventricular_min_qrs = ms_to_samples(120, fs)

    if qrslen > ventricular_min_qrs:
        return "V"

    if meta["p_pos"][channel] is not None:
        # qrs_start, qrs_end не могут быть None при штатном порядке вызова
        if qrslen < ventricular_min_qrs:
            return "S"
    else:
        if qrslen > ventricular_min_qrs:
            return "V"

    return "U"


def estimate_pq(meta):
    guess = []
    for i, p in enumerate(meta["p_pos"]):
        if p is None:
            continue
        q_pos =  meta["q_pos"][i]
        if q_pos is not None:
            guess.append(q_pos - p)
        else:
            r_pos = meta["r_pos"][i]
            if r_pos is not None:
                guess.append(r_pos - p)

    if guess:
        return np.mean(guess)


def estimate_qrslen(meta, fs, chan):
    """
    # Оценка длительности QRS-комплекса по зубцам
    :param meta:
    :param fs:
    :param chan:
    :return:
    """
    lb = int(meta["qrs_start"]*fs)
    rb = int(meta["qrs_end"]*fs)

    q_left = meta["q_pos"][chan]
    if q_left is not None:
        lb = q_left
    else:
        r_left = meta["r_start"][chan]
        if r_left is not None:
            lb = r_left

    s_right = meta["s_pos"][chan]
    if s_right is not None:
        rb = s_right
    else:
        r_right = meta["r_end"][chan]
        if r_right is not None:
            rb = r_right

    return rb - lb


def calculate_histogram(
        metadata,
        param_name,
        channel=None,
        nbins=10,
        censoring=False
):
    """
    Расчет гистограммы значений выбранного параметра в выбранном отведени
    :param metadata:
    :param param_name:
    :param channel:
    :param nbins:
    :param censoring: отбрасывание самых больших и самых маленьких значений
    :return:
    """

    if channel is None:
        param_val = [
            x[param_name] for x in metadata if x[param_name] is not None
        ]
    else:
        param_val = [
            x[param_name][channel] for x in metadata if x[param_name][channel] is not
            None
        ]

    # цензурирование - отбрасывем хвосты
    if censoring:
        q = np.percentile(param_val, [1,99])
        param_val = [x for x in param_val if q[0] < x < q[1]]

    hist, bine = np.histogram(param_val, bins=nbins)

    hdata = []

    for i, v in enumerate(hist):
        hdata.append(
            {
                "bin_left": bine[i],
                "bin_right": bine[i+1],
                "count": v,
                "percent": 100.0 * v / np.sum(hist)
            }
        )

    return hdata


def remove_outliers(metadata, key, dep_keys, delta):
    """
    Удаление далеко отстоящих точек
    :param x: массив значений
    :return: копия входного массива со значениями None для выбросов
    """

    xreal = [i for i in metadata[key] if i is not None]

    # нужно хотя бы три значения
    if len(xreal) < 3:
        return

    m = np.median(xreal)

    for i, v in enumerate(metadata[key]):
        if v is not None and abs(v - m) > delta:
            # print("{}[{}]: deviation {} samples".format(
            #     key, i, int(abs(v - m))))

            metadata[key][i] = None
            for k in dep_keys:
                metadata[k][i] = None


def detect_pvc(metadata, look=5, max_rr_s=1.5):
    """
    Поиск экстрасистол
    :param metadata:
    :param look: окрестность для построения модели
    :param max_rr_s: максимально допустимый RR-интервал
    :return: None. Изменяется поле flags в метаданных
    """

    # значащие RR-интервалы
    rr = []
    for i, qrs in enumerate(metadata):
        if not is_artifact(qrs):
            rri = qrs["RR"]
            if rri < max_rr_s:
                rr.append((i, qrs["qrs_center"], qrs["RR"]))

    # предварительно заполенный буфер RR-интервалов
    pack = rr[:look]

    for i, elem in enumerate(rr):

        if i >= len(rr)-1:
            break

        if i > look:
            pack.append(elem)

        if len(pack) > 2*look+1:
            pack.pop(0)

        slope, intercept, r_value, p_value, std_err1 = linregress(
            [x[1] for x in pack],
            [x[2] for x in pack]
        )

        # объединяем последний RR со следующим для проверки ЭС
        pack2 = pack[:]
        last = pack2.pop(-1)
        tst = rr[i+1]
        pack2.append((last[0], last[1], last[2] + tst[2]))

        slope, intercept, r_value, p_value, std_err2 = linregress(
            [x[1] for x in pack2],
            [x[2] for x in pack2]
        )

        if (std_err1 - std_err2) / std_err1 > 0.1:
            #print(tst[1], ["{:.2}".format(x[2]) for x in pack])
            metadata[tst[0]]["flags"] = "E"