#coding: utf-8

import sys
import os
import wfdb
import json
import numpy as np

if __name__ == "__main__":

    #filename = sys.argv[1]
    filename = '/Users/arseniy/SERDECH/data/PHYSIONET/I62'
    #write_dir = sys.argv[2]
    write_dir = '/Users/arseniy/SERDECH/data/converted'

    data, fields = wfdb.rdsamp(filename)
    rc = wfdb.rdrecord(
        filename,
        physical=False
    )

    #print(json.dumps(fields, indent=1))
    print(rc.adc_res)

    numch = data.shape[1]

    wfdb.wrsamp(
        os.path.basename(filename),
        fs=fields["fs"],
        units=fields["units"],
        sig_name=fields["sig_name"],
        p_signal=data,
        fmt=["16"] * numch,
        adc_gain=fields.get("adc_gain", [1.0]*numch),
        baseline=fields.get("baseline", [0]*numch),
        comments=fields["comments"],
        write_dir=write_dir
    )

    print(os.path.getsize(filename+".dat"))
    print(os.path.getsize(write_dir + "/" + os.path.basename(filename)
                          + ".dat"))
