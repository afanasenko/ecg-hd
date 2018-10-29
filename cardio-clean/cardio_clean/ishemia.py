# coding: utf-8

import numpy as np
from random import random, randint

from metadata import *
from config import config


ishemic_types = [
    "K1",  # нисходящая депрессия по критерию Kodama 1
    "K2",  # элевация по критерию Kodama 2
    "K3",  # элевация или депрессия по критерию Kodama 3
    "E1",  # нисходящая депрессия по критерию Ellestad 1
    "E2"   # косовосходящая депрессия по критерию Ellestad 2
]


def define_ishemia(metadata, chan, ishemic_type, start_idx, end_idx):

    start_time = metadata[start_idx]["qrs_start"]
    end_time = metadata[end_idx-1]["qrs_end"]

    hrs = [x["heartrate"] for x in metadata[start_idx:end_idx] if x[
        "heartrate"] is not None]
    sto = [abs(x["st_offset"][chan]) for x in metadata[start_idx:end_idx] if x[
        "st_offset"][chan] is not None]

    if hrs and sto:
        return {
            "type": ishemic_type,
            "channel": chan,
            "start": start_time,
            "end": end_time,
            "count": end_idx - start_idx + 1,
            "max_offset": np.max(sto),
            "heartrate": np.mean(hrs),
            "modified": False
        }


def define_ishemia_episodes(sig, header, metadata, **kwargs):
    """
        Поиск эпизодов ишемии по критериям Kodama и Ellestad в ST-сегментее
    :param sig:
    :param header:
    :param metadata:
    :param kwargs: См. секцию STT в config.yaml
    :return:
    """

    ishemia = []

    if not metadata:
        return ishemia

    numch = len(metadata[0]["r_pos"])
    last_codes = [None]*numch
    last_seq = [None]*numch
    fs = header["fs"]

    k2_dur = kwargs.get(
        "kodama_elev_dur",
        config.STT["kodama_elev_dur"]
    ) * fs
    e1_dur = kwargs.get(
        "ellestad_duration",
        config.STT["ellestad_duration"]
    ) * fs
    e2_dur = kwargs.get(
        "ellestad_duration",
        config.STT["ellestad_duration"]
    ) * fs
    min_len = kwargs.get(
        "min_episode",
        config.STT["min_episode"]
    )
    k1_duration = kwargs.get(
        "k1_duration",
        config.STT["k1_duration"]
    ) * fs

    for ch in range(numch):
        # перевод порогов из милливольт в значения сигнала
        bias = header["baseline"][ch]
        gain = header["adc_gain"][ch]
        k1_thresh = -kwargs.get(
            "kodama_depr_t",
            config.STT["kodama_depr_t"]
        ) * gain
        k2_thresh = kwargs.get(
            "kodama_elev_t",
            config.STT["kodama_elev_t"]
        ) * gain
        k3_thresh = kwargs.get(
            "kodama_relation",
            config.STT["kodama_relation"]
        ) * gain
        e1_thresh = -kwargs.get(
            "ellestad_depr_t1",
            config.STT["ellestad_depr_t1"]
        ) * gain
        e2_thresh = -kwargs.get(
            "ellestad_depr_t2",
            config.STT["ellestad_depr_t2"]
        ) * gain

        for i, meta in enumerate(metadata):

            if is_artifact(meta):
                continue

            c = ""

            stbeg = meta["st_start"][ch]
            stend = meta["st_end"][ch]
            jplus = meta["st_plus"][ch]
            jlev = meta["st_start_level"][ch]
            if jlev is not None:
                jlev *= gain
            stlev = meta["st_plus_level"][ch]
            if stlev is not None:
                stlev *= gain
            iso = meta["isolevel"][ch]
            if iso is not None:
                iso += bias
            slope = meta["st_slope"][ch]

            # критерии применимы только при положительном зубце T
            #if meta["t_height"][ch] is None or meta["t_height"][ch] < 0:
            #    continue

            # Kodama-1: горизонтальная либо нисходящая депрессия
            if stlev is not None and slope is not None\
                and stlev < k1_thresh and slope <= 0:
                    c = "K1"

            # Kodama-2: элевация не менее 80 мс от начала сегмента (J)
            elif jplus is not None and stbeg is not None and stend is not None\
                and jlev is not None and iso is not None and jlev > k2_thresh:

                elev = stbeg
                while sig[elev, ch] - iso > k2_thresh:
                    elev += 1
                    if elev >= stend:
                        break

                if elev - stbeg > k2_dur:
                    c = "K2"

            # Kodama-3: отношение смещения st к ЧСС больше порога (элевация или депр.)
            elif stlev is not None and meta["heartrate"] is not None\
                and abs(stlev) / meta["heartrate"] > k3_thresh:

                c = "K3"

            # Ellestad-1 косонисходящая депрессия
            elif jlev is not None and slope is not None:
                if jlev < e1_thresh and slope <= 0:
                    if jlev is not None and iso is not None:
                        depr = jlev

                        while sig[depr, ch] - iso < e1_thresh:
                            depr += 1
                            if depr >= stend:
                                break

                        if depr - stbeg > e1_dur:
                            c = "E1"

            # Ellestad-2 косовосходящая депрессия
            elif stlev is not None and slope is not None and stbeg is not \
                    None and stend is not None:
                if stlev < e2_thresh and slope > 0:
                    if stend - stbeg > e2_dur:
                        c = "E2"

            itype = last_codes[ch]

            if itype != c or i == len(metadata) - 1:

                if itype in ishemic_types and last_seq[ch] > min_len:

                    ep_start = i-last_seq[ch]
                    allow = True
                    if itype == "K1":
                        episode_duration = metadata[i]["qrs_end"] - \
                                           metadata[ep_start]["qrs_start"]
                        if episode_duration < k1_duration:
                            allow = False

                    if allow:
                        epi = define_ishemia(
                            metadata,
                            ch,
                            itype,
                            ep_start,
                            ep_start + last_seq[ch]
                        )
                        if epi is not None:
                            ishemia.append(epi)

                last_codes[ch] = c
                last_seq[ch] = 1
            else:
                last_seq[ch] += 1

    return ishemia


def mock_ishemia_episodes(metadata):
    """
    Псевдослучайная имитация ишемических жпизодов, для проверки интерфейса
    :param metadata:
    :return:
    """

    ishemia = []
    numch = len(metadata[0]["r_pos"])

    ish_prob = 0.5

    i = 0
    while i < len(metadata):

        chanid = randint(0, numch-1)
        typeid = randint(0, len(ishemic_types)-1)
        duration = randint(5, 100)
        end_cyc = min(len(metadata), i + duration)

        epi = define_ishemia(
            metadata,
            chanid,
            ishemic_types[typeid],
            i, end_cyc
        )

        if epi is not None and random() > ish_prob:
            ishemia.append(epi)

        i += duration

    return ishemia