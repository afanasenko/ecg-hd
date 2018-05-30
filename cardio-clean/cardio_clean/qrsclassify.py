# coding: utf-8

import numpy as np
from scipy.signal import lfilter

from sigbind import signal_channels


def correlate(sig1, sig2):
    """
        Вычисление статического коэффициента корреляции многомерных сигналов
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
    center = int(meta["r_wave_center"][0] * fs)
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
    # Номера классов - абстрактные, присваиваются исходя из количества
    # экземпляров
    for i, qcl in enumerate(
            sorted(
                qrs_classes,
                key=lambda x: len(x["samples"]),
                reverse=True
            )
    ):
        for s in qcl["samples"]:
            metadata[s]["qrs_class_id"] = i

            # помечаем ущербные классы как артефакты
            if len(qcl["samples"]) <= artifact_threshold:
                metadata[s]["artifact"] = True
            else:
                metadata[s]["artifact"] = False

        classdesc.append({
            "id": i,
            "average": qcl["accumulator"] / len(qcl["samples"]),
            "count": len(qcl["samples"])
        })

    # помечаем неклассифицированные комплексы как артефакты
    for m in metadata:
        if "qrs_class_id" not in m:
            m["qrs_class_id"] = None
            m["artifact"] = True

    return classdesc


def incremental_classifier(sig, hdr, metadata, classgen_t=0.9, include_data=0):
    """
        Однопроходный классификатор QRS-комплексов
    :param sig: сигнал
    :param hdr: заголовок сигнала
    :param metadata: метаданные с разметкой QRS-комплексов
    :param classgen_t: нижний порог коэффициента корреляции на создание нового класса
    :return: список классов
    """

    num_cyc = len(metadata)
    fs = hdr["fs"]

    first_qrs, first_c = extract_qrs(sig, fs, metadata[1])

    # инициализировать первый класс вторым комплексом, т.к. 1-й может быть
    # обрезан
    qrs_classes = [
        {
            "accumulator": first_qrs,
            "center": first_c,
            "samples": {1},
            "data": [(first_qrs, first_c)] if include_data else []
        }
    ]

    for i in range(num_cyc):

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
                        "samples": {i},
                        "data": [(new_qrs, new_c)] if include_data else []
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

                    if include_data and len(qcl["data"]) < include_data:
                        new_qrs, new_c = extract_qrs(sig, fs, metadata[i])
                        qcl["data"].append((new_qrs, new_c))

    # добавляем номера классов в метаданные
    # и формируем сокращенное описание каждого класса
    return finalize_classes(qrs_classes, metadata)
