# coding: utf-8

import numpy as np
from config import config

from metadata import estimate_pq, estimate_rr, estimate_pp, is_pvc
from rhythms import find_episodes


def is_rbbb(qrs):
    """
    Определяем наличие блокады ПНПГ
    :param qrs:
    :return: True/False
    """

    channel_seq = config.SIGNAL["default_channel_sequence"]

    # Нужны отведения V1 и V6, поэтому работает только на 12-канальном сигнале
    numch = len(qrs["r_pos"])
    if numch < 12:
        return False

    v1 = channel_seq.index("V1")

    # признак - отриц. T в V1
    tv1 = qrs["t_height"][v1]
    if tv1 is not None and tv1 < 0:
        return True

    # признак - раздвоенный R-зубец
    if qrs["r2_pos"][v1] is not None:
        return True

    # TODO: нужна доп. проверка на ширину QRS


def is_lbbb(qrs):
    """
    Определяем наличие блокады ЛНПГ
    :param qrs:
    :return: True/False
    """

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


def define_sablock(metadata):
    """
     Определяем SA-блокаду 2 и 3 степени
    :param metadata:
    :return: rythms
    """
    rythms = []
    total_cycles = len(metadata)
    prr = 0

    min_pause = 1.5
    max_pause_2 = 2.0

    for i, qrs in enumerate(metadata):
        if i:
            rr = qrs["RR"]
            if rr is None or is_pvc(qrs):
                continue

            if prr > 0 and rr > min_pause * prr and i < total_cycles-1:

                # это SA-блокада
                # в зависимости от продолжительности паузы назначаем тип
                blockade = "SAB_II" if rr <= max_pause_2*prr else "SAB_III"

                pause_start = qrs["qrs_center"]

                rythms.append({
                    "desc": blockade,
                    "start": pause_start,
                    "end": pause_start + rr,
                    "modified": False
                })

    return rythms


def define_avblock(metadata, fs, min_episode, **kwargs):
    """
     Определяем AV-блокаду
    :param metadata:
    :param fs: частота дискретизации
    :param kwargs: mobitz_pq, min_pause
    :return: rythms
    """

    rythms = []

    total_cycles = len(metadata)

    if not total_cycles:
        return rythms

    mobitz_pq = kwargs.get(
        "mobitz_pq",
        config.BLOCK["mobitz_pq"]
    ) * fs

    min_pause = kwargs.get(
        "min_pause",
        config.BLOCK["min_pause"]
    )

    # для оценки PQ используем окно wnd циклов
    wnd = 10
    avb_marks = [""] * total_cycles

    pqbuf = []
    ppbuf = []
    rrbuf = []

    avb_level = 0

    for ncycle, qrs in enumerate(metadata):

        if len(pqbuf) > wnd:
            pqbuf.pop(0)
        if len(ppbuf) > wnd:
            ppbuf.pop(0)
        if len(rrbuf) > wnd:
            rrbuf.pop(0)

        #
        pqi = estimate_pq(qrs)
        if pqi is not None:
            pqbuf.append(pqi)

        ppi = estimate_pp(metadata, ncycle)
        if ppi is not None:
            ppbuf.append(ppi)

        rri = estimate_rr(metadata, ncycle)
        if rri is not None:
            rrbuf.append(rri)

        if pqbuf:
            avg_pq = np.mean(pqbuf)

            if avg_pq > mobitz_pq:
                avb_level = 1

            # AV-блокада 2 степени Мобиц-2 - блуждающий PQ c выпадением цикла
            if np.std(pqbuf) > 0.1 * avg_pq:
                avb_level = 2

            if rrbuf and ppbuf:
                if np.mean(rrbuf) < 0.5 * np.mean(ppbuf):
                    avb_level = 3

        if avb_level > 1:
            if len(rrbuf) > 1 and rri > min_pause * rrbuf[-2]:

                pause_start = qrs["qrs_center"]

                rythms.append({
                    "desc": "AVB_II",
                    "start": pause_start,
                    "end": pause_start + rri,
                    "modified": False
                })
        else:
            if avb_level:
                # блокада 1-й и 3-й степени не сопровождается выпадением комплексов
                avb_marks[ncycle] = "AVB_" + "I" * avb_level

    # ищем эпизоды блокады 1-й степени
    rythms += find_episodes(avb_marks, min_episode, metadata)
    return rythms