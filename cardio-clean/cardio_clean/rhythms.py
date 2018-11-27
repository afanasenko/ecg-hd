# coding: utf-8

import numpy as np

from metadata import estimate_pq


rhythm_signatures = [
    ("sin_norm",  u"Нормальный синусовый ритм"),
    ("sin_tachy", u"Cинусовая тахикардия"),
    ("sin_brady", u"синусовая брадикардия"),
    ("sin_other",  u"Синусовая аритмия"),
    ("pacemaker_migration",  u"Миграция водителя ритма"),
    ("atrial",  u"Предсердный ритм"),
    ("ventricular",  u"Желудочковый ритм"),
    ("VT", u"Желудочковая тахикардия"),   # 100-400 bpm
    ("VTP",  u"Пароксизмальная желудочковая тахикардия"),
    ("av_parox", u"Пароксизмальная AB тахикардия"),
    ("a_parox", u"Пароксизмальная предсердная тахикардия"),
    ("PVC", u"Экстрасистолия"),  # можно убрать - разобрана на подтипы
    ("pause", u"Асистолия"),  # можно убрать - разобрана на подтипы
    ("LBBB", u"Блокада ЛНПГ"),
    ("RBBB", u"Блокада ПНПГ"),
    ("WPW", u"Синдром WPW"),
    ("AFL", u"Трепетание предсердий"),   # 200-400 bpm
    ("VFL", u"Трепетание желудочков"),   # 200-300 bpm
    ("AFB", u"Фибрилляция предсердий"),  # 350-700 bpm
    ("VFB", u"Фибрилляция желудочков"),   # 200-500 bpm
    ("PVC_SSE", u"Единичная наджелудочковая ранняя экстрасистолия"),
    ("PVC_SSI", u"Единичная наджелудочковая вставочная экстрасистолия"),
    ("PVC_SVE", u"Единичная желудочковая ранняя экстрасистолия"),
    ("PVC_SVI", u"Единичная желудочковая вставочная экстрасистолия"),
    ("PVC_CSE", u"Парная наджелудочковая ранняя экстрасистолия"),
    ("PVC_CSI", u"Парная наджелудочковая вставочная экстрасистолия"),
    ("PVC_CVE", u"Парная желудочковая ранняя экстрасистолия"),
    ("PVC_CVI", u"Парная желудочковая вставочная экстрасистолия"),
    ("PVC_GSE", u"Групповая наджелудочковая ранняя экстрасистолия"),
    ("PVC_GSI", u"Групповая наджелудочковая вставочная экстрасистолия"),
    ("PVC_GVE", u"Групповая желудочковая ранняя экстрасистолия"),
    ("PVC_GVI", u"Групповая желудочковая вставочная экстрасистолия"),
    ("SAB_I", u"СА-блокада 1 степени"),  # можно убрать - неразличима на ЭКГ
    ("SAB_II", u"СА-блокада 2 степени"),
    ("SAB_III", u"СА-блокада 3 степени"),
    ("AVB_I", u"АВ-блокада 1 степени"),
    ("AVB_II", u"АВ-блокада 2 степени"),
    ("AVB_III", u"АВ-блокада 3 степени")
]


def is_migration(metadata_block, pilot_chan):
    """
    Выявляем миграцию водителя ритма по смене полярности P-зубцов
    :param metadata_block:
    :param pilot_chan:
    :return: True/False
    """
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


def is_flutter(qrs):
    pilot_chan = 1 if len(qrs["r_pos"]) > 1 else 0
    # 0.3 - пороговое значение корреляции для обнаружения периодичности
    return qrs["flutter"][pilot_chan] > 0.3


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
                if count >= min_episode and last_r >= 0 and last_r != "":

                    start_time = metadata[i-count]["qrs_start"]
                    end_time = metadata[i-1]["qrs_end"]

                    desc = last_r if type(last_r) == str else rhythm_signatures[last_r][0]

                    rythms.append({
                        "desc": desc,
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