# coding: utf-8

"""Анализ работы каридиостимулятора
"""

import numpy as np

from cardio_clean.metadata import *
from cardio_clean.rhythms import find_episodes, rhythm_codes
from cardio_clean.config import config


def define_pacemaker_episodes(metadata, **kwargs):
    """Поиск эпизодов работы кардиостимулятора

    :param metadata:
    :param kwargs: min_episode
    :return:
    """

    min_episode = kwargs.get(
        "min_episode",
        config.PACEMAKER["min_episode"]
    )

    if not metadata:
        return []

    numch = len(metadata[0].get("pma",[]))
    if not numch:
        return []

    pm_marks = []
    total_cycles = len(metadata)

    for ch in range(numch):
        pm_marks.append(np.zeros(total_cycles, int) - 1)

    for ncycle, qrs in enumerate(metadata):
        for ch, pmach in enumerate(qrs["pma"]):
            if len(pmach) == 1:
                # проверяем положение спайка относительно P-зубца, если он есть
                ppos = safe_p_pos(qrs)
                if pmach[0] < ppos:
                    pm_marks[ch][ncycle] = rhythm_codes["PACED_A"]
                else:
                    pm_marks[ch][ncycle] = rhythm_codes["PACED_V"]
            elif len(pmach) == 2:
                pm_marks[ch][ncycle] = rhythm_codes["PACED_D"]

    # Ищем отведение, в котором стимулятор сильнее всего проявился
    selected_chan = None
    max_paced_cycles = 0
    for ch, rmarks in enumerate(pm_marks):
        paced = sum(rmarks >= 0)
        if paced > max_paced_cycles:
            max_paced_cycles = paced
            selected_chan = ch

    if selected_chan is not None:
        return find_episodes(
            pm_marks[selected_chan],
            min_episode=min_episode,
            metadata=metadata)
    else:
        return []