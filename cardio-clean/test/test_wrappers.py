# coding: utf-8

import os
from cardio_clean.cardioproc_api import *
from cardio_clean.metadata import samples_to_ms


def test_readwrite():
    filename_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "TestFromDcm.ecg")
    filename_out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "output_rw.ecg")

    with open(filename_in, "rb") as fi:
        hdr, data = read_buffer(fi)

    with open(filename_out, "wb") as fo:
        write_buffer(fo, hdr, data)

    assert os.stat(filename_in).st_size == os.stat(filename_out).st_size

    with open(filename_out, "rb") as fcheck:
        hdr2, data2 = read_buffer(fcheck)

        assert data.shape == data2.shape

    os.remove(filename_out)


def test_baseline():
    filename_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "TestFromDcm.ecg")
    filename_out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "output_bl.ecg")

    with open(filename_in, "rb") as fcheck:
        start_hdr, start_data = read_buffer(fcheck)

    with open(filename_in, "rb") as fi:
        with open(filename_out, "wb") as fo:
            blobapi_fix_baseline(inbuf=fi, outbuf=fo)

    assert os.stat(filename_in).st_size == os.stat(filename_out).st_size

    # Проверка размерности обработанного сигнала
    with open(filename_out, "rb") as fcheck:
        end_hdr, end_data = read_buffer(fcheck)

        assert start_data.shape == end_data.shape


    os.remove(filename_out)


def test_mains():
    filename_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "TestFromDcm.ecg")
    filename_out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "output_mains.ecg")

    with open(filename_in, "rb") as fi:
        with open(filename_out, "wb") as fo:
            blobapi_mains_correction(inbuf=fi, outbuf=fo)

    assert os.stat(filename_in).st_size == os.stat(filename_out).st_size

    with open(filename_out, "rb") as fcheck:
        hdr, data = read_buffer(fcheck)

        print("signal shape after mains correction: {}".format(data.shape))

    os.remove(filename_out)


def test_qrs():
    filename_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "TestFromDcm.ecg")

    with open(filename_in, "rb") as fi:
        meta = blobapi_detect_qrs(inbuf=fi)

    assert len(meta[0]) == 3628

    # Проверка на адекватность значения средней ЧСС
    avg_heartrate = np.median([x["heartrate"] for x in meta])

    assert 60.0 <= avg_heartrate <= 120


def test_parameters():
    filename_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "TestFromDcm.ecg")
    filename_out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "output_bl.ecg")

    with open(filename_in, "rb") as fi:
        with open(filename_out, "wb") as fo:
            blobapi_fix_baseline(inbuf=fi, outbuf=fo)

    with open(filename_out, "rb") as fi:
        meta = blobapi_detect_qrs(
            inbuf=fi,
            min_qrs_ms=20,
            channel=0,
            postprocessing=True
        )

        assert len(meta[0]) == 3628

        stdur = [nqrs["stt_params"]["duration"] for nqrs in meta[0] if nqrs["stt_params"]["duration"]]

        print("st-сегментов: {}\nСредняя длительность : {} мс".format(
            len(stdur),
            np.mean(stdur))
        )

        assert len(stdur) > 0


def disable_test_classify():

    filename_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "TestFromDcm.ecg")
    filename_out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "output_bl.ecg")

    print("Remove bias...")
    with open(filename_in, "rb") as fi:
        with open(filename_out, "wb") as fo:
            blobapi_fix_baseline(inbuf=fi, outbuf=fo)

    print("Find QRS...")
    with open(filename_out, "rb") as fi:
        meta = blobapi_detect_qrs(inbuf=fi)
        print("cycles found: {}".format(len(meta)))

    print("Classification...")
    with open(filename_out, "rb") as fi:
        classes = blobapi_classify_qrs(
            inbuf=fi,
            metadata=meta[0],
            classgen_threshold=0.8
        )

        print("classes found: {}".format(len(classes)))

        assert len(classes) == 1

        num_art = len([x for x in meta[0] if x["artifact"]])

        assert num_art == 0

    os.remove(filename_out)


if __name__ == "__main__":
    #test_readwrite()
    #test_baseline()
    #test_mains()
    test_parameters()
    #test_classify()