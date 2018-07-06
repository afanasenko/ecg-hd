import os
import datetime
from cardio_clean.cardioproc_api import *

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
        #print(meta)

    print("Classification...")
    with open(filename_out, "rb") as fi:
        print(meta[0])
        classes = blobapi_classify_qrs(
            inbuf=fi,
            metadata=meta,
            classgen_threshold=0.8
        )

        print("classes found: {}".format(len(classes)))
        print(classes)

        assert len(classes) == 1

        num_art = len([x for x in meta if x["artifact"]])

        assert num_art == 1

if __name__ == "__main__":
    print(datetime.datetime.now())
    test_classify()
    print(datetime.datetime.now())