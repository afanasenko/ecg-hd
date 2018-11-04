# coding: utf-8

import numpy as np
from random import randint
from config import config

from metadata import *


rhythm_signatures = [
    ("sin_norm",  u"Нормальный синусовый ритм"),
    ("sin_tachy", u"Cинусовая тахикардия"),
    ("sin_brady", u"синусовая брадикардия"),
    ("sin_other",  u"Синусовая аритмия"),
    ("pacemaker_migration",  u"Миграция водителя ритма"),
    ("atrial",  u"Предсердный ритм"),
    ("ventricular",  u"Желудочковый ритм"),
    ("v_parox",  u"Пароксизмальная желудочковая тахикардия"),
    ("av_parox", u"Пароксизмальная AB тахикардия"),
    ("a_parox", u"Пароксизмальная предсердная тахикардия"),
    ("PVC", u"Экстрасистолия"),
    ("pause", u"Асистолия"),
    ("LBBB", u"Блокада ЛНПГ"),
    ("RBBB", u"Блокада ПНПГ"),
    ("WPW", u"Синдром WPW"),
    ("AFL", u"Трепетание предсердий"),
    ("VFL", u"Трепетание желудочков"),
    ("AFB", u"Фибрилляция предсердий"),
    ("VFB", u"Фибрилляция желудочков")
]

# Цифровые коды ритмов
rhythm_names = {a:b[0] for a,b in enumerate(rhythm_signatures)}
rhythm_codes = {b[0]:a for a,b in enumerate(rhythm_signatures)}


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
        return u"Неизвестный ритм"


def is_migration(metadata_block, pilot_chan):
    pq = []
    pm = []
    for x in metadata_block:
        pqest = estimate_pq(x)
        if pqest is not None:
            pq.append(pqest)
        p = x["p_height"][pilot_chan]
        if p is not None:
            pm.append(p)

    if len(pq) > 0.5*len(metadata_block):
        if np.std(pq) > 0.5 * np.mean(pq):
            return True

    if len(pm) > 0.5*len(metadata_block):
        if np.std(pm) > 0.5 * np.mean(pm):
            return True

    return False


def define_pvc(metadata, rythms):

    for i, qrs in enumerate(metadata):
        if is_pvc(qrs) and not is_artifact(qrs):
            rythms.append({
                "desc": "PVC",
                "start": qrs["qrs_start"],
                "end": qrs["qrs_end"],
                "modified": False
            })


def define_pauses(metadata, rythms):
    total_cycles = len(metadata)
    for i, qrs in enumerate(metadata):
        if i:
            rr = qrs["RR"]
            if rr is None:
                continue

            prr = metadata[i-1]["RR"]
            if prr is not None and rr > 1.5 * prr and i < total_cycles-1:

                c_pause = (qrs["qrs_center"] + metadata[i+1]["qrs_center"]) / 2

                rythms.append({
                    "desc": "pause",
                    "start": c_pause - 0.1,
                    "end": c_pause + 0.1,
                    "modified": False
                })


def is_flutter(qrs):
    pilot_chan = 1 if len(qrs["r_pos"]) > 1 else 0
    #print(qrs["f_waves"][pilot_chan])
    return 2 < qrs["f_waves"][pilot_chan] < 10


def is_rbbb(qrs):

    channel_seq = config.SIGNAL["default_channel_sequence"]

    # Нужны отведения V1 и V6, поэтому работает только на 12-канальном сигнале
    numch = len(qrs["r_pos"])
    if numch < 12:
        return False

    v1 = channel_seq.index("V1")
    v6 = channel_seq.index("V6")

    # признак - отриц. T в V1
    tv1 = qrs["t_height"][v1]
    if tv1 is not None and tv1 < 0:
        return True

    # признак - раздвоенный R-зубец
    if qrs["r2_pos"][v1] is not None:
        return True

    # нужна доп. проверка на ширину QRS


def is_lbbb(qrs):

    channel_seq = config.SIGNAL["default_channel_sequence"]

    # Нужны отведения V1 и V6, поэтому работает только на 12-канальном сигнале
    numch = len(qrs["r_pos"])
    if numch < 12:
        return False

    v1 = channel_seq.index("V1")
    v6 = channel_seq.index("V6")

    # признак - отриц. T в V6
    tv6 = qrs["t_height"][v6]
    if tv6 is not None and tv6 < 0:
        return True

    # нужна доп. проверка на ширину QRS


def find_episodes(rythm_marks, min_episode, metadata):
    rythms = []
    last_r = None
    count = 0
    total_cycles = len(metadata)
    for i, r in enumerate(rythm_marks):
        if i:
            if r == last_r and i < total_cycles-1:
                count += 1
            else:
                if count >= min_episode and last_r >= 0:

                    start_time = metadata[i-count]["qrs_start"]
                    end_time = metadata[i-1]["qrs_end"]
                    rythms.append({
                        "desc": rhythm_signatures[last_r][0],
                        "start": start_time,
                        "end": end_time,
                        "modified": False
                    })
                last_r = r
                count = 1
        else:
            last_r = r
            count = 1
    return rythms


def define_rythm(metadata, **kwargs):
    """

    :param metadata:
    :param kwargs: См. секцию RHYTHM в config.yaml
    :return:
    """

    min_episode = kwargs.get(
        "min_episode",
        config.RHYTHM["min_episode"]
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
            rythm_marks[ncycle] = rhythm_codes["AFL"]
            continue

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

        num_sin = len([1 for x in ct if x == "N"])
        num_sup = len([1 for x in ct if x == "S"])
        num_ven = len([1 for x in ct if x == "V"])

        # постоянный RR-интервал
        if np.std(hr) < 0.1 * avg_hr:
            if 60.0 <= avg_hr <= 100.0:
                if num_sup > wnd:
                    rythm_marks[ncycle] = rhythm_codes["atrial"]
                else:
                    rythm_marks[ncycle] = rhythm_codes["sin_norm"]
            elif avg_hr < 60.0:
                if avg_hr <= 45.0 and num_ven > wnd:
                    rythm_marks[ncycle] = rhythm_codes["ventricular"]
                else:
                    rythm_marks[ncycle] = rhythm_codes["sin_brady"]
            else:
                rythm_marks[ncycle] = rhythm_codes["sin_tachy"]
        else:
            rythm_marks[ncycle] = rhythm_codes["sin_other"]

    rythms = find_episodes(rythm_marks, min_episode, metadata)
    rythms += find_episodes(syndrome_marks, min_episode, metadata)

    define_pvc(metadata, rythms)
    define_pauses(metadata, rythms)

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
