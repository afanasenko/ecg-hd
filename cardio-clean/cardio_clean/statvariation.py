# coding: utf-8

import numpy as np


class RecusiveStat():
    """utility class for recursive statistics calculation
    """

    def __init__(self):
        self.count = 0
        self.sum_val = 0.0
        self.sum_sq = 0.0

    def reset(self):
        self.count = 0
        self.sum_val = 0.0
        self.sum_sq = 0.0

    def size(self):
        return self.count

    def add(self, value):
        self.count += 1
        self.sum_val += value
        self.sum_sq += value*value

    def average(self):
        if self.count > 0:
            return self.sum_val / self.count
        else:
            return None

    def dispersion(self):
        count = self.count
        if count > 2:
            disper = self.sum_sq / count - self.sum_val*self.sum_val / (count*count)
            if disper >= 0:
                return disper
            else:
                return 0.0
        else:
            return None

    def std_dev(self):
        disp = self.dispersion()
        if disp is None:
            return None
        else:
            return np.sqrt(disp)


def rhythm_stats(metadata):
    """Расчет статистической вариабельности ритма

    :param metadata: список метаданных (только чтение). Требуемые ключи:
        - complex_type
        - RR

    :return: metrics: словарь (с фиксированным набором ключей), содержащий
    параметры вариабельностив виде имя-значение.
    Вычисляемые параметры: RRmean, SDNN, pNN50, SDANN, RMSSD, SDNNi.
    Все параметры имеют размерность [мс].
    При отсутствии достаточного числа нормальных
    qrs-комплексов все или часть параметров могут иметь значение None.
    """

    nn_stat = RecusiveStat()
    dnn_stat = RecusiveStat()
    ann_stat = RecusiveStat()
    five_stat = RecusiveStat()
    nni_stat = RecusiveStat()

    count50 = 0
    delta = 50      # порог на различие двух смежных интервалов
    last_nn = 0

    timestamp = metadata[0]["qrs_start"]

    for i, x in enumerate(metadata[:-1]):

        if x["qrs_start"] - timestamp > 300:
            # начало нового 5-минутного интервала
            if five_stat.size():
                five_av = five_stat.average()
                if five_av is not None:
                    ann_stat.add(five_av)
                five_std = five_stat.std_dev()
                if five_std is not None:
                    nni_stat.add(five_std)
                five_stat.reset()

            timestamp = x["qrs_start"]

        if x["complex_type"] == "N" and metadata[i+1]["complex_type"] == "N":
            # сразу переводим в мс
            rr = x["RR"]*1000

            if rr is not None:

                nn_stat.add(rr)
                five_stat.add(rr)

                if last_nn:
                    drr = rr - last_nn
                    if abs(drr) > delta:
                        count50 += 1

                    dnn_stat.add(drr)

                last_nn = rr

    rrm = nn_stat.average()
    rrsd = nn_stat.std_dev()
    pnn = 100.0 * count50 / nn_stat.size() if nn_stat.size() else None
    rmssd = dnn_stat.std_dev()
    sdann = ann_stat.std_dev()
    sdnni = nni_stat.average()

    return {
        "RRmean": rrm,
        "SDNN": rrsd,
        "pNN50": pnn,
        "SDANN": sdann,
        "RMSSD": rmssd,
        "SDNNi": sdnni
    }