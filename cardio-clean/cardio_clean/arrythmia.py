# coding: utf-8

import numpy as np
from random import randint

from metadata import *

rythm_signatures = [
    "sin_norm",  # - синусовый ритм: зубец Р есть во всех отведениях,
    # PII (+)
    "sin_tachy",
    "sin_brady",
    "sin_other",  # другая синусовая аритмия
    "pacemaker_migration",  # миграция водителя ритма
    "atrial",  # предсердный ритм: отрицательные PII, PIII, неизменный QRS
    "ventricular",  # желудочковый ритм: ЧСС<45, QRS>0,12 и не связан с P.
    "a_fib",  # фибрилляция предсердий
    "v_fib",  # трепетание желудочков
    "v_parox",  # пароксизмальная желудочковая тахикардия
    "av_parox",  # пароксизмальная наджелудочковая AB тахикардия
    "a_parox"  # пароксизмальная наджелудочковая предсердная тахикардия
]

# Цифровые коды ритмов
rythm_names = {a:b for a,b in enumerate(rythm_signatures)}
rythm_codes = {b:a for a,b in enumerate(rythm_signatures)}


def is_migration(metadata_block, pilot=1):
    pq = []
    pm = []
    for x in metadata_block:
        pqest = estimate_pq(x)
        if pqest is not None:
            pq.append(pqest)
        p = x["p_height"][pilot]
        if p is not None:
            pm.append(p)

    if len(pq) > 0.5*len(metadata_block):
        if np.std(pq) > 0.5 * np.mean(pq):
            return True

    if len(pm) > 0.5*len(metadata_block):
        if np.std(pm) > 0.5 * np.mean(pm):
            return True

    return False


def define_rythm(metadata):

    # для оценки средних показателей используем окно +- wnd циклов
    wnd = 10

    total_cycles = len(metadata)
    rythm_marks = np.zeros(total_cycles, int)

    for ncycle, qrs in enumerate(metadata):

        if ncycle < wnd:
            bnd = [0,2*wnd]
        elif ncycle > total_cycles-wnd:
            bnd = [total_cycles-2*wnd,total_cycles]
        else:
            bnd = [ncycle-wnd,ncycle+wnd]

        # обнаруживаем миграцию водителя
        if is_migration(metadata[bnd[0]:bnd[1]]):
            rythm_marks[ncycle] = rythm_codes["pacemaker_migration"]
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
                    rythm_marks[ncycle] = rythm_codes["atrial"]
                else:
                    rythm_marks[ncycle] = rythm_codes["sin_norm"]
            elif avg_hr < 60.0:
                if avg_hr <= 45.0 and num_ven > wnd:
                    rythm_marks[ncycle] = rythm_codes["ventricular"]
                else:
                    rythm_marks[ncycle] = rythm_codes["sin_brady"]
            else:
                rythm_marks[ncycle] = rythm_codes["sin_tachy"]
        else:
            rythm_marks[ncycle] = rythm_codes["sin_other"]

    rythms = []
    min_epi = wnd
    last_r = None
    count = 0
    for i, r in enumerate(rythm_marks):
        if i:
            if r == last_r and i < total_cycles-1:
                count += 1
            else:
                if count >= min_epi:

                    start_time = metadata[i-count]["qrs_start"]
                    end_time = metadata[i-1]["qrs_end"]
                    rythms.append({
                        "id": last_r,
                        "desc": rythm_signatures[last_r],
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


def mock_rythm_episodes(metadata):

    rythms = []

    i = 0
    while i < len(metadata):
        typeid = randint(0, len(rythm_codes)-1)
        duration = randint(20, 100)
        ending = min(len(metadata)-1, i + duration)

        start_time = metadata[i]["qrs_start"]
        end_time = metadata[ending]["qrs_end"]

        rythms.append({
            "id": typeid,
            "desc": rythm_signatures[typeid],
            "start": start_time,
            "end": end_time,
            "modified": False
        })

        i += duration

    return rythms
