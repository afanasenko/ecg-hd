# coding: utf-8

import numpy as np


def stt_analysis(x, fs, metadata):


    new_metadata = []
    valid_cycles = 0

    # guard
    min_st_samples = 3

    stt_params = {
        "st_start_offset": None,  # смещение начала ST от изолинии
        "st_end_offset": None,  # смещение конца ST от изолинии
        "st_slope": None,  # средний наклон ST
        "t_rising_slope": None,  # крутизна переднего фронта Т
        "t_falling_slope": None,  # крутизна заднего фронта Т
        "t_symmetry": None,  # симметричность Т (0 - идеальная симметрия)
        "t_sharpness": None,  # заостренность Т (чем больше, тем заостреннее)
        "t_amplitude": None  # амплитуда зубца T.
    }

    for ncycle, cycle_data in enumerate(metadata):

        params = cycle_data.copy()
        params["stt_params"] = stt_params.copy()

        if "waves" in cycle_data:
            st_start = cycle_data["waves"]["j"]["center"]
            st_end = cycle_data["waves"]["t"]["start"]

            if st_start is not None:
                if st_end is not None:
                    if st_end - st_start >= min_st_samples:

                        stt_params["st_start_offset"] = np.mean(
                            x[st_start:st_start+3]
                        )

                        stt_params["st_end_offset"] = np.mean(
                            x[st_end-3:st_end]
                        )

                        valid_cycles += 1

            new_metadata.append(params)

    return new_metadata
