#coding: utf-8

import numpy as np
from sigbind import fix_baseline, mains_filter
from qrsdetect import qrs_detection
from qrsclassify import incremental_classifier

# Числовые коды для поддерживаемых форматов
SAMPLE_TYPE_SHORT = 1  # 16-битный целочисленный со знаком

# когда в заголовке добавится поле SAMPLE_TYPE, этот флаг надо сбросить
legacy_headers=True


def read_buffer(buf):

    offset = 0
    count = 4 if legacy_headers else 5

    rawheader = np.fromfile(buf, dtype=np.uint32, count=count)

    signal_count = int(rawheader[0])
    signal_len = int(rawheader[1])
    sampling_frequency = int(rawheader[2])
    fp_bits = int(rawheader[3])

    sample_type = SAMPLE_TYPE_SHORT if legacy_headers else rawheader[4]

    if sample_type != SAMPLE_TYPE_SHORT:
        raise(Exception(("Unsupported sample type!")))

    offset += np.dtype(np.uint32).itemsize * count

    if fp_bits == 32:
        fptype = np.float32
    elif fp_bits == 64:
        fptype = np.float64
    else:
        raise

    adc_gain = np.fromfile(buf,
                             dtype=fptype,
                             count=signal_count)

    offset += np.dtype(fptype).itemsize * signal_count

    baseline = np.fromfile(buf,
                             dtype=np.uint32,
                             count=signal_count)

    offset += np.dtype(np.uint32).itemsize * signal_count

    data = np.fromfile(buf,
                         dtype=np.int16,
                         count=signal_count*signal_len)

    samples = np.reshape(np.array(data, float), (signal_len, signal_count))

    header = {
        "fs": sampling_frequency,
        "adc_gain": adc_gain,
        "baseline": baseline,
        "samples": signal_len,
        "channels": signal_count
    }

    return header, samples


def write_buffer(buf, header, samples):
    rawheader = np.array([
        samples.shape[1], # число каналов
        samples.shape[0], # число отсчетов
        header["fs"],
        header["adc_gain"].dtype.itemsize * 8,
    ], np.uint32)

    if not legacy_headers:
        header = np.append(header, SAMPLE_TYPE_SHORT)

    rawheader.tofile(buf)
    header["adc_gain"].tofile(buf)
    header["baseline"].tofile(buf)
    np.array(samples.flatten(), np.int16).tofile(buf)


def blobapi_fix_baseline(inbuf, outbuf, bias_window_ms=1500):
    """
       Обертка над функцией коррекции изолинии
    :param inbuf: входной буфер
    :param outbuf: модифициорованный буфер с заголовком и данными
    :param bias_window_ms: ширина окна для подаввления фона (мс)
    :return: None
    """

    header, indata = read_buffer(inbuf)
    outdata = fix_baseline(indata, header["fs"], bias_window_ms)
    write_buffer(outbuf, header, outdata)


def blobapi_mains_correction(
        inbuf,
        outbuf,
        attenuation=0.01,
        mains=50.0,
        aperture=512):
    """
       Обертка над функцией подавления сетевой помехи
    :param inbuf: входной буфер
    :param outbuf: модифициорованный буфер с заголовком и данными
    :param attenuation: коэффициент ослабления гармоник (0 - полное подавление)
    :param mains: частота сети
    :param aperture: апертура БПФ
    :return: None
    """

    header, indata = read_buffer(inbuf)
    outdata = mains_filter(
        indata,
        header["fs"],
        header["baseline"],
        mains=mains,
        attenuation=attenuation,
        aperture=aperture)
    write_buffer(outbuf, header, outdata)


def blobapi_detect_qrs(inbuf, min_qrs_ms=20):
    """
       Обнаружение QRS

        Предварительные условия:
        blobapi_detect_qrs желательно вызывать после функций подавления
        дрейфа изолинии blobapi_fix_baseline и подавления сетевой помехи
        blobapi_mains_correction (порядок вызова не имеет значения)

    :param inbuf: входной буфер (остается неизменным)
    :param min_qrs_ms: минимальная длительность QRS-комплекса
    :return: qrs_metadata (список найденных комплексов)
    """

    header, indata = read_buffer(inbuf)
    metadata, debugdata = qrs_detection(
        indata,
        fs=header["fs"],
        bias=header["baseline"],
        gain=header["adc_gain"],
        minqrs_ms=min_qrs_ms)

    return metadata


def blobapi_classify_qrs(inbuf, metadata, classgen_threshold=0.85):
    """
       Классификация QRS-комплексов

        Предварительные условия:
        blobapi_classify_qrs необходимо вызывать после того, как функцией
        blobapi_detect_qrs была выполнена разметка комплексов. Результаты
        классификации будут нболее точными, если сигнал был
        предварительно очищен от дрейфа изолинии и сетевой помехи.

        На вход функции вместе с сигналом подаются метаданные metadata,
        полученные из функции blobapi_detect_qrs. По результатам классификации
        в метаданные каждого комплекса добавляется поле qrs_class_id,
        содержащее порядковый номер класса для данного комплекса.

        В процессе работы функция разбивает множество обнаруженных
        QRS-комплексов на классы по подобию формы сигналов (учитываются все
        отведения). Число классов заранее неизвестно и регулируется параметром
        classgen_threshold. Чем больше значение порога, тем больше классов
        формируется в процессе разбиения. Нецелесообразно ставить порог
        выше 0.99 или ниже 0.25.

        Структура данных на выходе:
        qrs_classes = [
            {
                "id": порядковый номер класса (0 - самый часто встречающийся)
                "average": усредненные сигналы для класса (сохраняются все отведения)
                "count": число QRS-комплексов в классе
            },
            {
                "id",
                "average",
                "count"
            },
            ...
        ]

    :param inbuf: входной буфер (остается неизменным)
    :param metadata: метаданные ранее обнаруженных комплексов
    :param classgen_threshold: порог (от 0 до 1) на формирование новых
    классов (чем ниже порог, тем меньше получится классов)
    :return: список описаний найденных классов
    """

    header, indata = read_buffer(inbuf)
    qrs_classes = incremental_classifier(
        indata,
        header,
        metadata,
        classgen_t=classgen_threshold
    )

    return qrs_classes