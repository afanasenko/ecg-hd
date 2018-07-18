#coding: utf-8

# Цифровые коды ритмов
rcode = {
    "sinus": 0,  # - синусовый ритм: зубец Р есть во всех отведениях, PII (+)
#1 - атриовентрикулярный ритм
#2 - миграция водителя ритма
    "atrial": 3,  # - предсердный ритм: отрицательные PII, PIII, неизменный QRS
#4 - идиовентрикулярный (желудочковый внутрижелудочковый) ритм: QRS >
#160 (120) и не связан с P.
#5 - фибрилляция предсердий
#6 - трепетание желудочков
#7 - пароксизмальная желудочковая тахикардия
#8 - пароксизмальная наджелудочковая AB тахикардия
#9 - пароксизмальная наджелудочковая предсердная тахикардия
}


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
            qrs["rythm"] = rcode[guess]

