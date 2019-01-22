# coding: utf-8

from scipy.fftpack import fft
import numpy as np
from config import config
from metadata import is_artifact, is_pvc, safe_r_pos

from matplotlib import pyplot as plt


def rhythm_stats(metadata, fs):
    """Расчет статистической вариабельности ритма

    :param metadata: список метаданных (только чтение). Требуемые ключи:
        - complex_type
        - RR
    :param fs: частота дискретизации сигнала

    :return: metrics: словарь, содержащий параметры вариабельности
    в виде имя-значение. Вычисляемые параметры: RRmean, SDNN, pNN50, SDANN,
    RMSSD, SDNNi. При отсутствии достаточного числа нормальных
    qrs-комплексов все или часть параметров могут иметь значение None.
    Параметр SDANN содержит список, остальные - скалярные значения.
    """

    rrs = 0.0
    rrs2 = 0.0
    rrc = 0
    count50 = 0

    drrs = 0.0
    drrs2 = 0.0
    drrc = 0

    delta = 0.05*fs
    last_nn = 0

    timestamp = metadata[0]["qrs_start"]

    anns = 0.0
    anns2 = 0.0
    annc = 0

    for i, x in enumerate(metadata[:-1]):

        if x["qrs_start"] - timestamp > 300:

            if rrc:
                y = rrs/rrc
                annc += 1
                anns += y
                anns2 += y * y

            timestamp = x["qrs_start"]


        if x["complex_type"] == "N" and metadata[i+1]["complex_type"] == "N":
            rr = x["RR"]

            if rr is not None:

                rrc += 1
                rrs += rr
                rrs2 += rr*rr

                if last_nn:
                    drr = rr - last_nn
                    if abs(drr) > delta:
                        count50 += 1

                    drrc += 1
                    drrs += drr
                    drrs2 += drr*drr

                last_nn = rr

    if rrc:
        rrm = rrs / rrc
        rrsd = rrs2 / rrc - rrm * rrm
        pnn = 100.0 * count50 / rrc
    else:
        rrm = None
        rrsd = None
        pnn = None

    if drrc:
        rrm = drrs / drrc
        rmssd = drrs2 / drrc - rrm * rrm
    else:
        rmssd = None

    if annc:
        rrm = anns / annc
        sdann = anns2 / annc - rrm * rrm
    else:
        sdann = None

    return {
        "RRmean": rrm,
        "SDNN": rrsd,
        "pNN50": pnn,
        "SDANN": sdann,
        "RMSSD": rmssd,
        "SDNNi": None
    }