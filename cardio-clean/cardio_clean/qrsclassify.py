# coding: utf-8

"""Классификатор комплексов
"""

import numpy as np
from metadata import is_artifact, set_artifact, reset_artifact
from cardio_clean.config import config


def correlate(sig1, sig2):
    """ Вычисление статического коэффициента корреляции многомерных сигналов
    :param sig1:
    :param sig2:
    :return:
    """

    samples = sig1.shape[0]

    sx = 0.0
    sy = 0.0
    sxy = 0.0
    for i in range(samples):
        sxy += np.dot(sig1[i,:], sig2[i,:])
        sx += np.dot(sig1[i,:], sig1[i,:])
        sy += np.dot(sig2[i,:], sig2[i,:])

    return sxy / (np.sqrt(sx) * np.sqrt(sy))


def get_qrs_bounds(meta, fs):

    left = int(meta["qrs_start"] * fs)
    right = int(meta["qrs_end"] * fs)

    c = int(meta["qrs_center"] * fs)
    if c is None:
        center = int((left + right)/2)
    else:
        center = c
    return left, right, center


def qrs_correlation(sig, fs, meta0, meta1):

    l0, r0, c0 = get_qrs_bounds(meta0, fs)
    l1, r1, c1 = get_qrs_bounds(meta1, fs)
    leftw = min(c0-l0, c1-l1)
    rightw = min(r0-c0, r1-c1)

    qrs0 = sig[c0-leftw:c0+rightw, :]
    qrs1 = sig[c1-leftw:c1+rightw, :]

    cc = correlate(qrs0, qrs1)

    return float(cc)


def extract_qrs(sig, fs, meta):
    left, right, center = get_qrs_bounds(meta, fs)
    return sig[left:right, :], center - left


def qrs_reference_correlation(ref, sig, offset):

    smp = ref.shape[0]
    cc = correlate(
        ref,
        sig[offset:offset+smp, :]
    )

    return cc


def correlation_matrix(sig, hdr, metadata):

    num_cyc = len(metadata)
    corrmat = np.ones((num_cyc, num_cyc), float)
    fs = hdr["fs"]

    for i in range(num_cyc):
        for j in range(i+1, num_cyc):
            cc = qrs_correlation(sig, fs, metadata[i], metadata[j])
            corrmat[i,j] = cc

    print(np.min(corrmat))


def dominant_val(x):
    freq = {}
    for n in x:
        freq[n] = freq.get(n,0) + 1

    for k,v in sorted(freq.items(), key=lambda x: x[1], reverse=True):
        return k


def finalize_classes(qrs_classes, metadata):
    """
        Пост-обработка метаданных и вычисление усредненных комплексов
        помечает артефакты
    :param qrs_classes:
    :param metadata:
    :return:
    """

    artifact_threshold = 1
    classdesc = []
    known_complexes = set()

    # Номера классов - абстрактные, присваиваются исходя из количества
    # экземпляров
    for i, qcl in enumerate(
            sorted(
                qrs_classes,
                key=lambda x: len(x["samples"]),
                reverse=True
            )
    ):
        complex_types = []
        for s in qcl["samples"]:
            complex_types.append(metadata[s]["complex_type"])
            metadata[s]["qrs_class_id"] = i

            # помечаем ущербные классы как артефакты
            if len(qcl["samples"]) <= artifact_threshold:
                set_artifact(metadata[s])
            else:
                reset_artifact(metadata[s])
                known_complexes.add(s)

        if len(qcl["samples"]) > artifact_threshold:
            classdesc.append({
                "id": i,
                "average": qcl["accumulator"] / len(qcl["samples"]),
                "count": len(qcl["samples"]),
                "type": dominant_val(complex_types)
            })

    for i, qrs in enumerate(metadata):
        if i not in known_complexes:
            metadata[i]["qrs_class_id"] = None

    return classdesc


def incremental_classifier(sig, hdr, metadata, **kwargs):
    """ Однопроходный классификатор QRS-комплексов
    :param sig: сигнал (многоканальный)
    :param hdr: заголовок сигнала
    :param metadata: метаданные с разметкой QRS-комплексов
    :param kwargs: class_generation: нижний порог коэффициента корреляции на
    создание
    нового класса
    :return: список классов
    """

    classgen_t = kwargs.get(
        "class_generation",
        config.CLASSIFIER["class_generation"]
    )

    # число циклов
    num_cyc = len(metadata)

    if num_cyc < 2:
        print("Not enough data for classification")
        return {}

    fs = hdr["fs"]

    first_qrs, first_c = extract_qrs(sig, fs, metadata[1])

    # инициализируем первый класс вторым комплексом, т.к. 1-й может быть
    # обрезан
    qrs_classes = [
        {
            "accumulator": first_qrs,
            "center": first_c,
            "samples": {1}
        }
    ]

    for i in range(num_cyc):

        if is_artifact(metadata[i]):
            continue

        l1, r1, c1 = get_qrs_bounds(metadata[i], fs)
        cormat = np.zeros(len(qrs_classes), float)
        for c, qrsc in enumerate(qrs_classes):

            cand_offset = c1 - qrsc["center"]
            ref_samples = qrsc["accumulator"].shape[0]

            # пропускаем комплексы, присутствующие в сигнале не полностью
            if cand_offset < 0 or cand_offset+ref_samples > hdr["samples"]:
                continue

            cormat[c] = qrs_reference_correlation(
                qrsc["accumulator"] / len(qrsc["samples"]),
                sig,
                cand_offset
            )

        if np.any(cormat):
            max_cc = np.max(cormat)

            if max_cc < classgen_t:
                # создание нового класса
                new_qrs, new_c = extract_qrs(sig, fs, metadata[i])
                qrs_classes.append(
                    {
                        "accumulator": new_qrs,
                        "center": new_c,
                        "samples": {i}
                    }
                )
            else:
                # выбор существуюущего класса
                classid = np.argmax(cormat)
                qcl = qrs_classes[classid]

                # обновление усредненного цикла в выбранном классе
                left_acc = qcl["center"]
                right_acc = qcl["accumulator"].shape[0] - left_acc

                if c1 - left_acc >= 0 and c1 + right_acc <= sig.shape[0]:

                    qcl["accumulator"] += sig[c1-left_acc:c1+right_acc, :]
                    qcl["samples"].add(i)

    # добавляем номера классов в метаданные
    # и формируем сокращенное описание каждого класса
    return finalize_classes(qrs_classes, metadata)
