# coding: utf-8

import numpy as np

from metadata import is_artifact, set_pvc


def detect_pvc_episodes(metadata, fs):
    """
    Поиск экстрасистол
    :param metadata:
    :param look: окрестность для построения модели
    :param max_rr_s: максимально допустимый RR-интервал
    :return: None. Изменяется поле flags в метаданных
    """

    # сюда записываем все отметки всех циклов
    emarks = [""] * len(metadata)

    # порог обнаружения начала ЭС относительно нормального интервала RR
    scep = 0.65

    # общая длительность вставочной ЭС относительно нормального интервала RR
    ins = 1.2

    # значащие RR-интервалы
    rr = [x["RR"] for x in metadata if not is_artifact(x)]
    # предварительная оценка нормального RR
    rr_e = np.median(rr)*fs

    i = 0
    pvc_state = 0
    sdur = 0
    egroup = []

    while i < len(metadata):
        qrs = metadata[i]
        if not is_artifact(qrs):

            rr_smp = qrs["RR"] * fs

            if pvc_state == 0:
                if rr_smp < scep * rr_e and i < len(metadata) - 1:
                    pvc_state = 1

            if pvc_state == 1:
                sdur = rr_smp
                pvc_state = 2
            else:
                if pvc_state == 2:
                    sdur += rr_smp
                    egroup.append(i)

                    k = 1.01
                    if sdur > k*rr_e or rr_smp >= scep*rr_e:
                        # разбираемся с накопившимися экстрасистолами
                        # ранняя или вставочная
                        atr1 = "I" if sdur < ins * rr_e else "E"
                        for ncycle in egroup:
                            # желудочковая или наджелудочковая
                            atr2 = metadata[ncycle]["complex_type"]
                            if atr2 != "V":
                                atr2 = "S"

                            # одиночная/парная/групповая
                            atr3 = "G"
                            if len(egroup) == 1:
                                atr3 = "S"
                            elif len(egroup) == 2:
                                atr3 = "C"

                            emarks[ncycle] = "PVC_" + atr3 + atr2 + atr1
                            set_pvc(metadata[ncycle])

                        pvc_state = 0
                        sdur = 0
                        egroup = []

        i += 1

    #print(".".join(emarks))
    return emarks