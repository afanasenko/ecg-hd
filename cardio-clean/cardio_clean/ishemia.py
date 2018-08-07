# coding: utf-8

import numpy as np
from random import random, randint

from metadata import *


ishemic_types = [
    "K1",  # нисходящая депрессия по критерию Kodama 1
    "K2",  # элевация по критерию Kodama 2
    "K3",  # элевация или депрессия по критерию Kodama 3
    "E1",  # нисходящая депрессия по критерию Ellestad 1
    "E2"   # косовосходящая депрессия по критерию Ellestad 2
]


def check_criteria(sig, chan, meta, gain, fs):

    # перевод порогов из милливольт в значения сигнала
    k1_thresh = 0.1 * gain
    k2_thresh = 0.1 * gain
    k2_dur = 0.08 * fs
    k3_thresh = 0.0014 * gain
    e1_thresh = -0.1 * gain
    e1_dur = 0.08 * fs
    e2_thresh = -0.2 * gain
    e2_dur = 0.08 * fs

    stbeg = meta["st_start"][chan]
    stend = meta["st_end"][chan]
    jplus = meta["st_plus"][chan]
    jlev = meta["st_start_level"][chan]
    stlev = meta["st_plus_level"][chan]
    iso = meta["isolevel"][chan]

    # Kodama-1: горизонтальная либо нисходящая депрессия
    if stlev is not None and meta["st_slope"] is not None:
        if stlev < -k1_thresh and meta["st_slope"] <= 0:
            return "K1"

    # Kodama-2: элевация не менее 80 мс от начала сегмента (J)
    if jplus is not None and stbeg is not None and stend is not None:
        if jlev is not None and iso is not None:

            if jlev > k2_thresh:
                elev = stbeg

                while sig[elev, chan] - iso > k2_thresh:
                    elev += 1
                    if elev >= stend:
                        break

                if elev - stbeg > k2_dur:
                    return "K2"

    # Kodama-3: отношение смещения st к ЧСС больше порога (элевация или депр.)
    if stlev is not None and meta["heartrate"] is not None:
        if abs(stlev) / meta["heartrate"] > k3_thresh:
            return "K3"

    # Ellestad-1
    if jlev is not None and meta["st_slope"] is not None:
        if jlev < e1_thresh and meta["st_slope"] <= 0:
            if jlev is not None and iso is not None:
                    depr = jlev

                    while sig[depr, chan] - iso < e1_thresh:
                        depr += 1
                        if depr >= stend:
                            break

                    if depr - stbeg > e1_dur:
                        return "E1"
    # Ellestad-2
    if stlev is not None and meta["st_slope"] is not None and stbeg is not \
            None and stend is not None:
        if stlev < e2_thresh and meta["st_slope"] > 0:
            if stend - stbeg > e2_dur:
                return "E2"

    return ""


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


def define_ishemia_episodes(sig, header, metadata):

    ishemia = []
    numch = len(metadata[0]["r_pos"])
    last_codes = [None]*numch
    last_seq = [None]*numch
    min_len = 5
    k1_duration = 60.0 * header["fs"]

    for ch in range(numch):

        for i, meta in enumerate(metadata):

            if meta["artifact"]:
                continue

            c = check_criteria(
                sig,
                ch,
                meta,
                header["adc_gain"][ch],
                header["fs"]
            )

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