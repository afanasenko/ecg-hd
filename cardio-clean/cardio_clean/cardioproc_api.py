#coding: utf-8

import numpy as np
from sigbind import fix_baseline, mains_filter
from qrsdetect import qrs_detection

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

    samples = np.reshape(np.array(data, float), (-1,signal_count))

    header = {
        "fs": sampling_frequency,
        "adc_gain": adc_gain,
        "baseline": baseline,
    }

    return header, samples


def write_buffer(buf, header, samples):
    rawheader = np.array([
        samples.shape[1],
        samples.shape[0],
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
        mains=mains,
        attenuation=attenuation,
        aperture=aperture)
    write_buffer(outbuf, header, outdata)


def blobapi_detect_qrs(inbuf, min_qrs_ms=20):
    """
       Обертка над функцией обнаружения QRS
    :param inbuf: входной буфер (остается неизменным)
    :param min_qrs_ms: минимальная длительность QRS-комплекса
    :return: qrs_metadata (список найденных комплексов)
    """

    header, indata = read_buffer(inbuf)
    metadata, debugdata = qrs_detection(
        indata,
        header["fs"],
        minqrs_ms=min_qrs_ms)

    return metadata
