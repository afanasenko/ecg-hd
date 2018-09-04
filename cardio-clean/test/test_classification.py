# coding: utf-8

import os
import datetime
from cardio_clean.cardioproc_api import *
from cardio_clean.metadata import is_artifact


def test_classify():

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
        meta = blobapi_detect_qrs(inbuf=fi, min_qrs_ms=20, postprocessing=True)
        print("cycles found: {}".format(len(meta)))
        print(meta[0])

    print("Classification...")
    with open(filename_out, "rb") as fi:

        classes = blobapi_classify_qrs(
            inbuf=fi,
            metadata=meta,
            classgen_threshold=0.8
        )

        print("classes found: {}".format(len(classes)))
        print(classes)

        assert len(classes) == 1

        artifacts = [i for i,x in enumerate(meta) if is_artifact(x)]
        print(artifacts)
        num_art = len(artifacts)

        assert num_art > 1

        num_classified = len(
            [x for x in meta if x["qrs_class_id"] is not None]
        )

        classids = set(
            [x["qrs_class_id"] for x in meta]
        )

        print(num_classified)
        print(classids)

    with open(filename_out, "rb") as fi:
        ry = blobapi_classify_rythms(
            inbuf=fi,
            metadata=meta,
        )

        assert len(ry) > 1

    with open(filename_out, "rb") as fi:
        ish = blobapi_find_ishemia(
            inbuf=fi,
            metadata=meta,
        )

        assert len(ish) == 0

    os.remove(filename_out)


if __name__ == "__main__":
    print(datetime.datetime.now())
    test_classify()
    print(datetime.datetime.now())