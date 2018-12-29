# coding: utf-8

"""
.. module:: util
   :platform: Unix, Windows
   :synopsis: Вспомогательные функции для чтения/записи сигналов ЭКГ
Модуль содержит вспомогательные функции для чтения/записи сигналов ЭКГ
"""

import wfdb
import numpy as np


# Числовые коды для поддерживаемых форматов
SAMPLE_TYPE_SHORT = 1  # 16-битный целочисленный со знаком

# когда в заголовке добавится поле SAMPLE_TYPE, этот флаг надо сбросить
legacy_headers=True


def read_buffer(buf):
    """Чтение сигнала из открытого двоичного файла

    :param buf: Файловый объект
    :type buf: fileobj
    :returns: dict - заголовок файла
    :returns: numpy.array - Массив отсчетов
    """
    offset = 0
    count = 4 if legacy_headers else 5

    rawheader = np.fromfile(buf, dtype=np.uint32, count=count)

    signal_count = int(rawheader[0])
    signal_len = int(rawheader[1])
    sampling_frequency = int(rawheader[2])
    fp_bits = int(rawheader[3])

    sample_type = SAMPLE_TYPE_SHORT if legacy_headers else rawheader[4]

    if sample_type != SAMPLE_TYPE_SHORT:
        raise(Exception("Unsupported sample type!"))

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


def signal_channels(signal):
    """Перебор каналов многоканального ЭКС

    :param signal: np.array массив с размерностью (samples, channels)
    :return: объект-генератор возвращает № канала и массив отсчетов
    """
    if len(signal.shape) == 1:
        yield 0, np.array(signal, np.float32)
    else:
        numch = signal.shape[1]
        for chan in range(numch):
            yield chan, np.array(signal[:,chan], np.float32)


def ecgread(filename):
    """Чтение сигнала из файла (поддерживаются форматы ecg и MIT-BIH)

    :param filename: имя файла (в случае MIT-BIH передается без расширения)
    :type filename: str
    :return:
    """
    if filename.endswith(".ecg"):
        with open(filename, "rb") as fi:
            hdr, data = read_buffer(fi)
            print(hdr)
            return data, hdr

    else:
        data, fields = wfdb.rdsamp(filename)
        # rdsamp возвращает сигнал без смещения в физических единицах
        numch = data.shape[1]
        hdr = {
            "fs": fields["fs"],
            "adc_gain": np.array([1.0]*numch),
            "baseline": np.array([0.0]*numch),
            "samples": data.shape[0],
            "channels": data.shape[1]
        }

        print(fields["sig_name"])


        return data, hdr

common_signal_names = ['I', 'II', 'III', 'AVR', 'AVL', 'AVF', 'V1', 'V2',
                       'V3', 'V4', 'V5', 'V6']