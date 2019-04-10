# coding: utf-8

import numpy as np
from random import randint
from config import config

from rhythms import *
from pvcdetect import detect_pvc_episodes
from blockades import *
from pmdetect import define_pacemaker_episodes


def arrythmia_name(s):
    """
    Получение человекочитаемого названия для заданного ритма
    :param s: краткое название ритма
    :return: длинное название ритма на русском
    """
    try:
        n = rhythm_signatures.index(s)
        return rhythm_signatures[n][1]
    except:
        return u"Неизвестный ритм " + unicode(s)


def define_rythm(metadata, **kwargs):
    """

    :param metadata:
    :param kwargs: min_episode, fs
    :return:
    """

    min_episode = kwargs.get(
        "min_episode",
        config.RHYTHM["min_episode"]
    )

    fs = kwargs.get(
        "fs",
        config.SIGNAL["default_fs"]
    )

    if not metadata:
        return []

    # для оценки средних показателей используем окно +- wnd циклов
    wnd = 10
    pilot_chan = 1 if len(metadata[0]["r_pos"]) > 1 else 0

    total_cycles = len(metadata)
    rythm_marks = np.zeros(total_cycles, int)
    # синдромы могут перекрываться с основными ритмами, поэтоу храним отдельно
    # FIXME: хранить все вместе
    syndrome_marks = np.zeros(total_cycles, int) - 1
    flutter_marks = np.zeros(total_cycles, int) - 1

    for ncycle, qrs in enumerate(metadata):

        if is_lbbb(qrs):
            syndrome_marks[ncycle] = rhythm_codes["LBBB"]
        elif is_rbbb(qrs):
            syndrome_marks[ncycle] = rhythm_codes["RBBB"]

        if ncycle < wnd:
            bnd = [0,2*wnd]
        elif ncycle > total_cycles-wnd:
            bnd = [total_cycles-2*wnd,total_cycles]
        else:
            bnd = [ncycle-wnd,ncycle+wnd]

        # обнаруживаем миграцию водителя
        if is_migration(metadata[bnd[0]:bnd[1]], pilot_chan):
            rythm_marks[ncycle] = rhythm_codes["pacemaker_migration"]
            continue

        # обнаруживаем трепетания во II-м отведении
        if is_flutter(qrs):
            flutter_marks[ncycle] = rhythm_codes["AFL"]

        # ЧСС
        hr = np.array(
            [x["heartrate"] for x in metadata[bnd[0]:bnd[1]] if x[
            "heartrate"] is not None]
        )

        avg_hr = np.mean(hr)

        # Доминирующий тип комплексов
        ct = np.array(
            [x["complex_type"] for x in metadata[bnd[0]:bnd[1]] if x[
            "complex_type"] != "U"]
        )

        dominant = "N"
        num_sin = len([1 for x in ct if x == "N"])
        num_sup = len([1 for x in ct if x == "S"])
        num_ven = len([1 for x in ct if x == "V"])
        if num_sup > wnd:
            dominant = "S"
        elif num_ven > wnd:
            dominant = "V"

        # постоянный RR-интервал
        r_mark = "sin_norm"
        if np.std(hr) < 0.1 * avg_hr:
            if 60.0 <= avg_hr <= 100.0:
                if dominant == "S":
                    r_mark = "atrial"
                elif dominant == "V":
                    # для желудочкового ритма это уже тахикардия
                    r_mark = "VT"

            elif avg_hr < 60.0:
                if dominant == "V":
                    r_mark = "ventricular"
                else:
                    r_mark = "sin_brady"
            else:
                r_mark = "sin_tachy"
        else:
            if 200 < avg_hr < 300 and dominant == "V":
                r_mark = "VFL"
            else:
                r_mark = "sin_other"

        rythm_marks[ncycle] = rhythm_codes[r_mark]

    rythms = find_episodes(rythm_marks, min_episode, metadata)
    rythms += find_episodes(syndrome_marks, min_episode, metadata)
    rythms += find_episodes(flutter_marks, min_episode, metadata)

    # сначала выделяем все ЭС
    pvc_marks = detect_pvc_episodes(metadata, fs)
    # потом размечаем их как обычные эпизоды, но с мин. длительностью 1
    rythms += find_episodes(pvc_marks, min_episode=1, metadata=metadata)
    # в последнюю очередь ищем паузы
    rythms += define_sablock(metadata)
    rythms += define_avblock(metadata, fs, min_episode)

    # и кардиостимулятор
    rythms += define_pacemaker_episodes(metadata)

    return rythms


def mock_rythm_episodes(metadata):

    rythms = []

    i = 0
    while i < len(metadata):
        typeid = randint(0, len(rhythm_codes)-1)
        duration = randint(20, 100)
        ending = min(len(metadata)-1, i + duration)

        start_time = metadata[i]["qrs_start"]
        end_time = metadata[ending]["qrs_end"]

        rythms.append({
            "desc": rhythm_signatures[typeid],
            "start": start_time,
            "end": end_time,
            "modified": False
        })

        i += duration

    return rythms
