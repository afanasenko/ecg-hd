# coding: utf-8

"""
Метаданные подразделяются на первичные и вторичные.
Первичные являются результатом автоматической сегментации и могут быть
скорректированы вручную. Для расчета вторичных данных используются как
первичные, так и сам сигнал. После ручного редактирования необходимо
пересчитывать вторичные данные, не проводя повторной сегментации сигнала.

# ######################################
Первичные метаданные включают:

Q, R, S, R* : центр зубца (метки времени)
P, Т: начало, центр, конец (метки времени)

qrsType: вид комплекса. Обнаруживаемые разновидности qrs - R, qR, Rs, qs, qRs.
Артефакты кодируются значением qrsType==None

# ######################################

Вторичные данные по зубцам
Q, R, S, R* (амплитуда)
P: амплитуда
T: амплитуда, крутизна переднего фронта, крутизна заднего фронта, симметрия, острота

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


def metadata_new():
    return {
        "waves": {
            "p": {"start": None, "end": None, "center": None, "height": None},
            "q": {"start": None, "end": None, "center": None, "height": None},
            "r": {"start": None, "end": None, "center": None, "height": None},
            "s": {"start": None, "end": None, "center": None, "height": None},
            "t": {"start": None, "end": None, "center": None, "height": None}
        },
        "qrsType": None,
        "RR": None,
        "heartrate": None,
        "isolevel": None,
        "ST": {
            "start": None,
            "stplus": None,
            "end": None,
            "start_level": None,
            "stplus_level": None,
            "end_level": None,
            "offset": None,
            "duration": None,
            "slope": None
        }
    }


def samples_to_ms(smp, fs):
    return smp * 1000.0 / fs


def ms_to_samples(ms, fs):
    return int(ms * fs / 1000.0)


def metadata_postprocessing(metadata, sig, fs, **kwargs):
    """
        Расчет вторичных параметров сигнала в одном отведении
    :param metadata:
    :param sig:
    :param fs:
    :param kwargs:
    :return:
    """

    j_offset_ms = kwargs.get("j_offset", 60)
    jplus_offset_ms = kwargs.get("jplus_offset", 80)

    for ncycle, cycledata in enumerate(metadata):

        # ######################################
        # точки J и J+
        # ставится со смещением от R-зубца
        rc = cycledata["waves"]["r"]["center"]
        if rc is not None:
            j_point = rc + ms_to_samples(j_offset_ms, fs)
            if j_point > len(sig) - 1:
                j_point = None
                jplus_point = None
            else:
                jplus_point = j_point + ms_to_samples(jplus_offset_ms, fs)
                if jplus_point > len(sig) - 1:
                    jplus_point = None

                elif cycledata["waves"]["t"]["center"] is not None:
                    tstart = cycledata["waves"]["t"]["center"]
                    jplus_point = min(j_point, tstart)

        else:
            j_point = None
            jplus_point = None

        # ######################################
        # запись высоты зубцов
        bias = cycledata.get("isolevel", 0.0)
        for wave in cycledata["waves"]:
            pos = cycledata["waves"][wave].get("center", None)
            if pos is None:
                cycledata["waves"][wave]["height"] = None
            else:
                cycledata["waves"][wave]["height"] = sig[pos] - bias

        # ######################################
        # RR
        if ncycle:
            cur_r = cycledata["waves"]["r"]["center"]
            prev_r = metadata[ncycle-1]["waves"]["r"]["center"]
            if all((cur_r, prev_r)):
                rr = samples_to_ms(cur_r - prev_r, fs)
                cycledata["RR"] = rr
                cycledata["heartrate"] = 60000.0 / rr

        # ######################################
        # ST
        st_start = j_point
        st_plus = jplus_point
        st_end = cycledata["waves"]["t"]["start"]

        cycledata["ST"]["start"] = st_start
        cycledata["ST"]["stplus"] = st_plus
        cycledata["ST"]["end"] = st_end

        if st_start is None:
            cycledata["ST"]["start_level"] = None
        else:
            cycledata["ST"]["start_level"] = sig[st_start] - bias

        if st_plus is None:
            cycledata["ST"]["stplus_level"] = None
        else:
            cycledata["ST"]["stplus_level"] = sig[st_plus] - bias

        if st_end is None:
            cycledata["ST"]["end_level"] = None
        else:
            cycledata["ST"]["end_level"] = sig[st_end] - bias

        if all((st_start, st_end)):
            dur = samples_to_ms(st_end - st_start, fs)
            if dur > kwargs.get("min_st_duration", 80):
                cycledata["ST"]["duration"] = dur
                cycledata["ST"]["offset"] = np.mean(sig[st_start:st_end]) - bias
                cycledata["ST"]["slope"] = (cycledata["ST"]["end_level"] -
                                        cycledata["ST"]["start_level"]) / cycledata["ST"]["duration"]
            else:
                cycledata["ST"]["duration"] = None


