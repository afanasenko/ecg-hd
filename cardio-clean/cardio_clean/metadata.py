# coding: utf-8


def metadata_new():
    return {
        "waves": {
            "p": {"start": None, "end": None, "center": None, "height": None},
            "q": {"start": None, "end": None, "center": None, "height": None},
            "r": {"start": None, "end": None, "center": None, "height": None},
            "s": {"start": None, "end": None, "center": None, "height": None},
            "j": {"center": None, "height": None},
            "j+": {"center": None, "height": None},
            "t": {"start": None, "end": None, "center": None, "height": None}
        },
        "qrsType": None,
        "RR": None,
        "ST": {
            "start_level": None,
            "end_level": None,
            "duration": None,
            "slope": None
        }
    }


def samples_to_ms(smp, fs):
    return smp * 1000.0 / fs


def metadata_postprocessing(metadata, sig, fs, **kwargs):
    """

    :param metadata:
    :param sig:
    :param fs:
    :param kwargs:
    :return:
    """

    j_offset_ms = kwargs.get("j_offset", 60)
    jplus_offset_ms = kwargs.get("jplus_offset", 80)

    for ncycle, cycledata in enumerate(metadata):

        # ######################################
        # точки J и J+
        # ставится со смещением от R-зубца
        rc = cycledata["waves"]["r"]["center"]
        if rc is not None:
            j_point = rc + int(fs*j_offset_ms/1000.0)
            if j_point > len(sig) - 1:
                j_point = None
                jplus_point = None
            else:
                jplus_point = j_point + int(fs * jplus_offset_ms / 1000.0)
                if jplus_point > len(sig) - 1:
                    jplus_point = None

        else:
            j_point = None
            jplus_point = None

        cycledata["waves"]["j"] = {"center": j_point}
        cycledata["waves"]["j+"] = {"center": jplus_point}

        # ######################################
        # запись высоты зубцов
        for wave in cycledata["waves"]:
            pos = cycledata["waves"][wave].get("center", None)
            if pos is not None:
                cycledata["waves"][wave]["height"] = sig[pos]

        # ######################################
        # RR
        if ncycle:
            cur_r = cycledata["waves"]["r"]["center"]
            prev_r = metadata[ncycle-1]["waves"]["r"]["center"]
            if all((cur_r, prev_r)):
                cycledata["RR"] = samples_to_ms(cur_r - prev_r, fs)
