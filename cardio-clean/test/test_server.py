#coding: utf-8

def test_qrs_result(input_file, kod_elev, kod_depr, elles_elev, elles_depr,
                 min_episide):

    header, indata = write_blob_to_memory(input_file)

    baseline_out = fix_baseline(indata, header["fs"], 1500)

    mains_out = mains_filter(baseline_out,
                             header["fs"],
                             header["baseline"],
                             mains=50,
                             attenuation=0.05,
                             aperture=512)

    metadata = qrs_detection(mains_out,
                             fs=header["fs"],
                             minqrs_ms=20)[0]

    find_points(mains_out,
                fs=header["fs"],
                metadata=metadata,
                bias=header["baseline"],
                gain=header["adc_gain"],
                debug=False)

    metadata_postprocessing(metadata,
                            mains_out,
                            header)

    try:
        qrs_classes = incremental_classifier(mains_out,
                                            header,
                                            metadata,
                                                classgen_t=0.998)
    except Exception as exception_obj:
        qrs_classes = []

    try:
        rythms = define_rythm(metadata, fs=header["fs"])
    except Exception as exception_obj:
        rythms = []

    try:
        ishemia = define_ishemia_episodes(mains_out,
                                          header,
                                          metadata,
                                          kodama_elev_t=kod_elev,
                                          kodama_depr_t=kod_depr,
                                          ellestad_depr_t1=elles_elev,
                                          ellestad_depr_t2=elles_depr,
                                          min_episode=min_episide)
    except Exception as exception_obj:
        ishemia = []

    classes_count = len(qrs_classes)
    for i in range(classes_count):
        qrs_classes[i]["average"] = np.reshape(qrs_classes[i]["average"],
                                               (1, np.product(qrs_classes[i]["average"].shape)))[0].tolist()

    try:
        s_var = rhythm_spectrum(metadata)
    except Exception as exception_obj:
        s_var = []

    try:
        stat_var = rhythm_stats(metadata)
    except Exception as exception_obj:
        stat_var = []



    return metadata, qrs_classes, rythms, ishemia, s_var, stat_var
