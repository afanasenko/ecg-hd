# coding: utf-8

from sigbind import signal_channels

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
        "qrs_class_id": None,
        "artifact": True,  # bool
        "qrsType": None,  # string

        # отдельные зубцы
        "p_start": [None]*num_channels,  # int array
        "p_end": [None]*num_channels,  # int array
        "p_pos": [None]*num_channels,  # int array
        "p_height": [None]*num_channels,  # float array
        "q_pos": [None]*num_channels,  # int array
        "q_height": [None]*num_channels,  # float array
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
        "st_slope": [None]*num_channels  # float array
    }


def samples_to_ms(smp, fs):
    return smp * 1000.0 / fs


def ms_to_samples(ms, fs):
    return int(ms * fs / 1000.0)


def level_from_pos(d, chan, pos_key, val_key, sig, bias):
    pos = d[pos_key][chan]
    if pos is None:
        d[val_key][chan] = None
    else:
        d[val_key][chan] = sig[pos] - bias


def metadata_postprocessing(metadata, sig, fs, **kwargs):
    """
    Расчет вторичных параметров сигнала в одном отведении

    Поскольку источник входных метаданных неизвестен, необходимо
    перезаписать значения всех зависимых ключей.
    :param metadata:
    :param sig:
    :param fs:
    :param kwargs: константы j_offset, jplus_offset_ms, min_st_ms
    :return: None (результатом являются измененные значения в metadata)
    """

    j_offset_ms = kwargs.get("j_offset", 60)
    jplus_offset_ms = kwargs.get("jplus_offset", 80)

    numch = sig.shape[1] if sig.ndim == 2 else 1

    # ритм оценивается всегда по второму отведению
    heartbeat_channel = 1 if numch > 1 else 0

    for ncycle, cycledata in enumerate(metadata):

        for chan, x in signal_channels(sig):

            # ######################################
            # точки J и J+
            # ставится со смещением от R-зубца
            rc = cycledata["r_pos"][chan]
            if rc is not None:
                j_point = rc + ms_to_samples(j_offset_ms, fs)
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

            level_from_pos(cycledata, chan, "p_pos", "p_height", x, bias)
            level_from_pos(cycledata, chan, "q_pos", "q_height", x, bias)
            level_from_pos(cycledata, chan, "r_pos", "r_height", x, bias)
            level_from_pos(cycledata, chan, "s_pos", "s_height", x, bias)
            level_from_pos(cycledata, chan, "t_pos", "t_height", x, bias)
            level_from_pos(cycledata, chan, "st_start",
                           "st_start_level", x, bias)
            level_from_pos(cycledata, chan, "st_plus", "st_plus_level", x,
                           bias)
            level_from_pos(cycledata, chan, "st_end", "st_end_level", x, bias)

            # ######################################
            # ST (продолжение)
            if all((j_point, st_end)):
                dur = samples_to_ms(st_end - j_point, fs)
                if dur > kwargs.get("min_st_ms", 80):
                    cycledata["st_duration"][chan] = dur

                    cycledata["st_offset"][chan] = np.mean(
                        sig[j_point:st_end]) - bias

                    cycledata["st_slope"][chan] = \
                        (cycledata["st_end_level"][chan] -
                         cycledata["st_start_level"][chan]) / \
                        cycledata["st_duration"][chan]
                else:
                    cycledata["st_duration"][chan] = None

            # ######################################
            # RR
            if chan == heartbeat_channel:
                if ncycle:
                    neighbour = metadata[ncycle-1]["r_pos"][chan]
                else:
                    neighbour = metadata[ncycle+1]["r_pos"][chan]

                if rc is not None and neighbour is not None:
                    rr = samples_to_ms(abs(rc - neighbour), fs)
                    cycledata["RR"] = rr
                    cycledata["heartrate"] = 60000.0 / rr
