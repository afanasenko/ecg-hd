# coding: utf-8

from cardio_clean.metadata import *
import numpy as np
from matplotlib import pyplot as plt
import json


def prepare_linplan(npoints):
    """ Расчет матрицы плана эксперимента для линейной регрессии

    :param npoints: число точек
    :return: матрица [npoints X 2]
    """
    p = np.ones((npoints,2))
    for i in xrange(npoints):
        p[i,0] = i

    pp = np.dot(np.transpose(p), p)

    return np.dot(np.linalg.pinv(pp), np.transpose(p))


def calculate_to_ts(rr_buf, pvc_pos, regres_window=4):

    sum_rr_left = rr_buf[pvc_pos-2] + rr_buf[pvc_pos-3]
    sum_rr_right = rr_buf[pvc_pos+2] + rr_buf[pvc_pos+1]
    to = 100.0 * (sum_rr_right - sum_rr_left) / sum_rr_left

    ts = 0

    pln = prepare_linplan(regres_window)

    for i in range(pvc_pos+1, len(rr_buf)-regres_window):

        cf = np.dot(pln, rr_buf[i:i+regres_window])
        slope = cf[0]

        #slope, intercept, r_value, p_value, std_err1 = linregress(
        #    [n for n in range(i, i + regres_window + 1)],
        #    rr_buf[i:i+regres_window+1]
        #)
        if slope > ts:
            ts = slope

    return to, ts


def turbulence_analyse(metadata, **kwargs):
    """ Анализ турбулентности ритма

    :param metadata: метаданные
    :param kwargs: beats_before, beats_after
    :return: turbulence_data - данные об эпизодах, отобранных для анализа
             mean_turbulence_data - усредненные показатели турбулентности
    """

    lh = int(kwargs.get(
        "beats_before",
        config.HRT["beats_before"]
    ))
    rh = int(kwargs.get(
        "beats_after",
        config.HRT["beats_after"]
    ))
    tswnd = int(kwargs.get(
        "regression_window",
        config.HRT["regression_window"]
    ))

    # сразу в мс
    rrbuf = [1000.0*qrs["RR"] if qrs["RR"] is not None else 0 for qrs in \
            metadata]
    pvc = [is_pvc(qrs) for qrs in metadata]

    turbulence_data = []

    # для тренда
    trend = None
    to_buf = []
    ts_buf = []

    for i, qrs in enumerate(metadata):
        if lh <= i < len(metadata)-rh:
            if pvc[i]:
                #  для анализа нужны фрагменты без экстрасистол
                # TODO: добавить проверку на тахи/брадикардию
                if any(pvc[i-lh:i]) or any(pvc[i+1:i+rh+1]):
                    continue

                if all(rrbuf[i-lh:i+rh]):
                    to, ts = calculate_to_ts(rrbuf, i, regres_window=tswnd)
                    turbulence_data.append(
                        {
                            "qrs_index": i,
                            "start_index": -lh + 1,
                            "TO": to,
                            "TS": ts,
                            "curve": rrbuf[i-lh:i+rh]
                        }
                    )
                    # для усреднения
                    to_buf.append(to)
                    ts_buf.append(ts)
                    if trend is None:
                        trend = np.array(rrbuf[i-lh:i+rh])
                    else:
                        trend += rrbuf[i-lh:i+rh]

    if turbulence_data:
        mean_turbulence_data = {
            "start_index": -lh+1,
            "TO": np.mean(to_buf),
            "TS": np.mean(ts_buf),
            "curve": trend / len(turbulence_data)
        }
    else:
        mean_turbulence_data = {
            "start_index": 0,
            "TO": None,
            "TS": None,
            "curve": np.array([])
        }

    return turbulence_data, mean_turbulence_data


if __name__ == "__main__":

    #metaname = "/Users/arseniy/SERDECH/data/PHYSIONET/I24.json"
    metaname = "/Users/arseniy//Downloads/Test20191007.ecg.json"

    print("Load...")
    meta = json.load(open(metaname))
    print("{} cycles".format(len(meta)))
    print("Start...")
    turb_data, trend_data = turbulence_analyse(meta)

    print([x["qrs_index"] for x in turb_data])
    print("TO={}, TS={}".format(trend_data["TO"], trend_data["TS"]))

    x0 = trend_data["start_index"]
    x1 = x0 + len(trend_data["curve"])

    xx = [i for i in range(x0, x1)]

    for turb in turb_data:
        plt.plot(xx, turb["curve"], 'b', alpha=0.3)

    #print(json.dumps(turb_data[0:2], indent=1))

    if trend_data["curve"] is not None:
        plt.plot(xx, trend_data["curve"], "r")

    plt.xticks(xx)

    plt.show()