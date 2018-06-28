# coding: utf-8

import numpy as np
from metadata import ms_to_samples, samples_to_ms


def stt_analysis(x, fs, metadata):

    new_metadata = []

    # guard value
    min_st_samples = ms_to_samples(80.0, fs)

    stt_params = {
        "start": None,  # метка времени начала начала ST
        "end": None,  # метка времени конца ST
        "start_level": None,  # смещение начала ST от изолинии
        "end_level": None,  # смещение конца ST от изолинии
        "offset": None,  # общее смещение ST от изолинии
        "duration": None,  # длительность ST в мс
        "slope": None,  # средний наклон ST
        #"t_rising_slope": None,  # крутизна переднего фронта Т
        #"t_falling_slope": None,  # крутизна заднего фронта Т
        #"t_symmetry": None,  # симметричность Т (0 - идеальная симметрия)
        #"t_sharpness": None,  # заостренность Т (чем больше, тем заостреннее)
        #"t_amplitude": None  # амплитуда зубца T.
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

                        bias = cycle_data.get("isolevel", 0.0)

                        stt_params["start"] = st_start
                        stt_params["end"] = st_end

                        stt_params["duration"] = samples_to_ms(
                            st_end - st_start, fs)

                        stt_params["start_level"] = x[st_start] - bias
                        stt_params["end_level"] = x[st_end] - bias

                        stt_params["offset"] = np.mean(
                            x[st_start:st_end]
                        ) - bias

                        stt_params["slope"] = \
                            (stt_params["end_level"] - stt_params[
                                "start_level"]) / stt_params["duration"]

            new_metadata.append(params)

    return new_metadata
