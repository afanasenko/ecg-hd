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

        assert len(classes) == 2

        artifacts = [i for i,x in enumerate(meta) if is_artifact(x)]
        print(artifacts)
        num_art = len(artifacts)

        assert num_art == 1

        num_classified = len(
            [x for x in meta if x["qrs_class_id"] is not None]
        )

        classids = set(
            [x["qrs_class_id"] for x in meta]
        )

        print("ClassIDs: {}, Number of classified complexes: {}".format(
            classids, num_classified))


    with open(filename_out, "rb") as fi:
        ry = blobapi_classify_rythms(
            inbuf=fi,
            metadata=meta,
        )

        assert len(ry) > 1

        pauses = [x for x in ry if x["desc"] == "pause"]
        print("Asystolia episodes:")
        print(pauses)

    with open(filename_out, "rb") as fi:
        ish = blobapi_find_ishemia(
            inbuf=fi,
            metadata=meta,
        )

        assert len(ish) == 0

    os.remove(filename_out)


def test_lowlevel():
    filename_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "testFindPoint.ecg")

    with open(filename_in, "rb") as f:
        header, indata = read_buffer(f)

    baseline_out = fix_baseline(indata, header["fs"], 1500)

    mains_out = mains_filter(baseline_out,
                             header["fs"],
                             header["baseline"],
                             mains=50,
                             attenuation=0.05,
                             aperture=512)

    print('start qrs_detection')
    metadata = qrs_detection(mains_out,
                             fs=header["fs"],
                             minqrs_ms=20)[0]

    print('start find_points')
    find_points(mains_out,
                fs=header["fs"],
                metadata=metadata,
                bias=header["baseline"],
                debug=False)

    print('start metadata_postprocessing')
    metadata_postprocessing(metadata,
                            mains_out,
                            header)

    print('start incremental_classifier')
    qrs_classes = incremental_classifier(mains_out,
                                         header,
                                         metadata,
                                         classgen_t=0.998)

    print('start define_rythm')
    rythms = define_rythm(metadata)

    print('start define_ishemia_episodes')
    ishemia = define_ishemia_episodes(mains_out,
                                      header,
                                      metadata,
                                      kodama_elev_t=0.1,
                                      kodama_depr_t=0.1,
                                      ellestad_depr_t1=0.2,
                                      ellestad_depr_t2=0.1,
                                      min_episode=5)

    print('start reshape classes average')
    classes_count = len(qrs_classes)
    for i in range(classes_count):
        qrs_classes[i]["average"] = np.reshape(qrs_classes[i]["average"],
                                               (1, np.product(qrs_classes[i][
                                                                  "average"].shape)))[
            0].tolist()


if __name__ == "__main__":
    print(datetime.datetime.now())
    test_lowlevel()
    print(datetime.datetime.now())