#coding: utf-8

from util import *
from sigbind import fix_baseline, mains_filter
from qrsdetect import qrs_detection
from qrsclassify import incremental_classifier
from wavdetect import find_points
from metadata import metadata_postprocessing, calculate_histogram
from arrythmia import define_rythm
from ishemia import define_ishemia_episodes

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


def blobapi_detect_qrs(
        inbuf,
        min_qrs_ms=20,
        postprocessing=True
    ):
    """
       Обнаружение QRS и последующая сегментация

        Предварительные условия:
        blobapi_detect_qrs желательно вызывать после функций подавления
        дрейфа изолинии blobapi_fix_baseline и подавления сетевой помехи
        blobapi_mains_correction (порядок вызова не имеет значения)

    :param inbuf: входной буфер (остается неизменным)
    :param min_qrs_ms: минимальная длительность QRS-комплекса
    :param postprocessing: расчет вторичных параметров (ритм, ST и др.)
    :return: [meta0, meta1, ..., metaN], где metaN - список найденных
    комплексов в N-ом канале
    """

    header, indata = read_buffer(inbuf)
    metadata = qrs_detection(
        indata,
        fs=header["fs"],
        minqrs_ms=min_qrs_ms
    )[0]

    # поиск характерных точек производится во всех отведениях
    find_points(
        indata,
        fs=header["fs"],
        metadata=metadata,
        bias=header["baseline"],
        debug=False
    )
    if postprocessing:
        metadata_postprocessing(
            metadata,
            indata,
            header
        )
    return metadata


def blobapi_postprocessing_qrs(
        inbuf,
        metadata
    ):
    """
    Расчет вторичных параметров на основе ранее
    сегментированного сигнала

        Предварительные условия:
        сигнал д.б. предварительно сегментирован функцией blobapi_detect_qrs

    :param inbuf:
    :param metadata: список словарей с данными сегментации каждого комплекса
    :return: None (изменяется содержимое metadata)
    """

    header, indata = read_buffer(inbuf)

    metadata_postprocessing(
        metadata,
        indata,
        header
    )


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


def blobapi_classify_rythms(
        inbuf,
        metadata
):
    """
        Выделение эпизодов нарушения ритма
    :param inbuf:
    :param metadata: метаданные с результатами сегментации
    :return: список эпизодов нарушения ритма
    [
     {
            "id": (int) цифровое обозначение аритмии,
            "desc": (string) текстовое описание
            "start": (float) начало эпизода в секундах от начала записи
            "end": (float) конец эпизода в секундах от начала записи
            "modified": (bool) флаг ручного редактирования
     },
     ...
    ]
    Мнимальная длительность эпизода принята равной 10 циклам для всех
    видов ритма. Поэтому фрагмент записи может не принадлежать ни одному из
    эпизодов, если в нем были обнаружены короткие эпизобы различных ритмов.
    """

    return define_rythm(metadata)


def blobapi_find_ishemia(
            inbuf,
            metadata
    ):
    """
        Выделение эпизодов, содержащих признаки ишемии в ST-сегменте
        Предварительные условия: выполнена классификация (blobapi_classify_qrs)
    :param inbuf:
    :param metadata: метаданные с результатами сегментации
    :return: список эпизодов ишемии
    [
     {
            "type": (string) тип эпизода ишемии, K1|K2|K3|E1|E2
            "channel": (int) номер отведения
            "start": (float) начало эпизода в секундах от начала записи
            "end": (float) конец эпизода в секундах от начала записи
            "count": (int) число комплексов данного типа
            "max_offset": (float) максимальное смещение ST в эпизоде
            "heartrate": (float) средняя ЧСС в эпизоде
            "modified": (bool) флаг ручного редактирования
     },
     ...
    ]
    """

    header, indata = read_buffer(inbuf)
    return define_ishemia_episodes(indata, header, metadata)


def blobapi_histogram_qt(metadata, channel=1):
    """
        Расчет гистограммы длин интервалов QT
    :param metadata: метаданные с результатами сегментации
    :param channel: интересующий номер отведения
    :return: данные для построения гистограммы в виде списка интервальных
    значений
    [
     {
        "bin_left": левая граница интервала (в секундах)
        "bin_right": правая граница интервала (в секундах)
        "count": число попаданий в интервал
        "percent": процент попаданий в интервал от общего числа (0 - 100)
     },
     ...
    ]
    Элементы списка уже упорядочены по положению интервала
    """
    return calculate_histogram(metadata, "qt_duration", channel)


def blobapi_histogram_qtc(metadata, channel=1):
    """
        Расчет гистограммы длин интервалов корригированного QT
    :param metadata: метаданные с результатами сегментации
    :param channel: интересующий номер отведения
    :return: см. описание для blobapi_histogram_qt
    """
    return calculate_histogram(metadata, "qtc_duration", channel)