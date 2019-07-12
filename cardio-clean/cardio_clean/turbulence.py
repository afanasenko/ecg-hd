# coding: utf-8

from cardio_clean.metadata import *
from scipy.stats import linregress
from matplotlib import pyplot as plt
import json


def calculate_to_ts(rr_buf, pvc_pos, regres_window=4):

    delta_rr_left = rr_buf[pvc_pos-2] - rr_buf[pvc_pos-3]
    delta_rr_right = rr_buf[pvc_pos+2] - rr_buf[pvc_pos+1]
    to = 100.0 * (delta_rr_right - delta_rr_left) / delta_rr_right

    ts = 0
    for i in range(pvc_pos+1, len(rr_buf)-regres_window):

        slope, intercept, r_value, p_value, std_err1 = linregress(
            [n for n in range(i, i + regres_window + 1)],
            rr_buf[i:i+regres_window+1]
        )
        if slope > ts:
            ts = slope

    return to, ts


def turbulence_analyse(metadata):
    """ Анализ турбулентности ритма

    :param metadata: метаданные
    :return: turbulence_data - данные об эпизодах, отобранных для анализа
             mean_turbulence_data - усредненные показатели турбулентности
    """

    # сразу в мс
    rrbuf = [1000.0*qrs["RR"] for qrs in metadata]
    pvc = [is_pvc(qrs) for qrs in metadata]

    lh = 4
    rh = 15
    tswnd = 4

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
            "TO": np.mean(to_buf),
            "TS": np.mean(ts_buf),
            "curve": trend / len(turbulence_data)
        }
    else:
        mean_turbulence_data = {
            "TO": None,
            "TS": None,
            "curve": None
        }

    return turbulence_data, mean_turbulence_data


if __name__ == "__main__":

    meta = json.load(open("/Users/arseniy/SERDECH/data/PHYSIONET/I13.json"))
    turb_data, trend_data = turbulence_analyse(meta)

    print([x["qrs_index"] for x in turb_data])
    print("TO={}, TS={}".format(trend_data["TO"], trend_data["TS"]))

    for turb in turb_data:
        plt.plot(turb["curve"], 'b', alpha=0.3)

    print(json.dumps(turb_data[0:2], indent=1))

    plt.plot(trend_data["curve"], "r")

    plt.show()