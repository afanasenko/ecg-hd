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