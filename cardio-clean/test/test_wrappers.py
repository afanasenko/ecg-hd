# coding: utf-8

import os
import json
from cardio_clean.cardioproc_api import *
from cardio_clean.metadata import is_artifact


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

    assert len(meta) == 3628

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

    with open(filename_in, "rb") as fi:
        meta = blobapi_detect_qrs(
            inbuf=fi,
            min_qrs_ms=20,
            postprocessing=True
        )

        filename_metadump = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "TestFromDcm.json"
        )

        with open(filename_metadump, "w") as fmeta:
            json.dump(meta, fmeta, indent=1, sort_keys=True)

        assert len(meta) == 3628

        # Проверка на наличие ЧСС для всех компелксов кроме артефактов
        no_rr = {c for c, x in enumerate(meta) if x["heartrate"] is None}
        arti = {c for c, x in enumerate(meta) if is_artifact(x)}

        assert arti == no_rr

        stdur = [x["st_duration"][0] for x in meta if x[
            "st_duration"][0] is not None]

        print(
            "{} циклов, {} st-сегментов\nСредняя длительность: {} мс".format(
                len(meta),
                len(stdur),
                np.mean(stdur)
            )
        )

        assert len(stdur) > 0

        h = blobapi_histogram_qt(metadata=meta, histogram_bins=11)
        h = blobapi_histogram_qtc(metadata=meta, histogram_bins=11)

    with open(filename_in, "rb") as fi:
        ish = blobapi_find_ishemia(inbuf=fi, metadata=meta,
                                   kodama_elev_dur=0.09)

    os.remove(filename_out)


if __name__ == "__main__":
    test_parameters()

