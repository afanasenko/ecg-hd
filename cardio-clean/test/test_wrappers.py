import os
from cardio_clean.wrappers import *


def test_readwrite():
    filename_in = "TestFromDcm.ecg"
    filename_out = "output_rw.ecg"

    with open(filename_in, "rb") as fi:
        hdr, data = read_buffer(fi)

        print(hdr["fs"])
        print(data.shape)

    with open(filename_out, "wb") as fo:
        write_buffer(fo, hdr, data)

    assert os.stat(filename_in).st_size == os.stat(filename_out).st_size


def test_baseline():
    filename_in = "TestFromDcm.ecg"
    filename_out = "output_bl.ecg"

    with open(filename_in, "rb") as fi:
        with open(filename_out, "wb") as fo:
            blobapi_fix_baseline(inbuf=fi, outbuf=fo)

    assert os.stat(filename_in).st_size == os.stat(filename_out).st_size


def test_mains():
    filename_in = "TestFromDcm.ecg"
    filename_out = "output_mains.ecg"

    with open(filename_in, "rb") as fi:
        with open(filename_out, "wb") as fo:
            blobapi_mains_correction(inbuf=fi, outbuf=fo)

    assert os.stat(filename_in).st_size == os.stat(filename_out).st_size

if __name__ == "__main__":
    test_readwrite()
    test_baseline()
    test_mains()
