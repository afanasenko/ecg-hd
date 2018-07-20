#coding: utf-8

from random import randint

rythm_signatures = [
    "sin_norm",  # - синусовый ритм: зубец Р есть во всех отведениях,
    # PII (+)
    "sin_tachy",
    "sin_brady",
    "sin_other",
# атриовентрикулярный ритм
# миграция водителя ритма
    "p_migration",
    "atrial",  # - предсердный ритм: отрицательные PII, PIII, неизменный QRS
    "ventricular",  #4 - идиовентрикулярный (желудочковый внутрижелудочковый)
    # ритм: QRS >
#160 (120) и не связан с P.
# фибрилляция предсердий
    "a_fib",
# трепетание желудочков
    "v_fib",
# - пароксизмальная желудочковая тахикардия
    "v_parox",
# - пароксизмальная наджелудочковая AB тахикардия
    "av_parox",
# - пароксизмальная наджелудочковая предсердная тахикардия
    "a_parox"
]

# Цифровые коды ритмов
rythm_codes = {b:a for a,b in zip(rythm_signatures, range(len(rythm_signatures)))}

def is_sinus_rythm(meta):
    """
    Проверка на синусовый ритм по положению и полярности P-зубцов
    :param meta:
    :param chan:
    :return:
    """

    std2 = 1
    std3 = 2

    # R-зубцы присутствуют во всех отведениях
    if all(meta["r_pos"]):
        for i, p in enumerate(meta["p_pos"]):
            # P-зубцы присутствуют во всех отведениях
            if p is None:
                return None
            # P-зубцы идут перед R-зубцами
            if meta["r_pos"][i] < p:
                return None

        # P-зубец во II-м отведении  положительный
        if meta["p_height"][std2] > 0:
            return "sinus"
        # P-зубец во II-м и III-м отв. отрицательный
        elif meta["p_height"][std3] < 0:
            return "atrial"


def find_rythm(sig, metadata):

    for ncycle, qrs in metadata:
        guess = is_sinus_rythm(qrs)
        if guess is not None:
            qrs["rythm"] = rythm_codes[guess]


def mock_rythm_episodes(metadata):

    rythms = []

    i = 0
    while i < len(metadata):
        type = randint(0, len(rythm_codes)-1)
        duration = randint(20, 100)
        ending = min(len(metadata)-1, i + duration)

        start_time = metadata[i]["qrs_start"]
        end_time = metadata[ending]["qrs_end"]

        rythms.append({
            "id": type,
            "вуыс": rythm_codes[type],
            "start": start_time,
            "end": end_time,
            "modified": False
        })

    return rythms
