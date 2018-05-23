# coding: utf-8

import io
import os
import sys
import numpy as np
from matplotlib import pyplot as plt, rc
import wfdb
from dtureader import dturead

if __name__ == "__main__":
    
    if len(sys.argv) > 1:
        fn = sys.argv[1]
        
        d, fs1, n = dturead(fn)
        signal = d
        
        n_sig = 6
        wfdb.wrsamp("name", fs = 4000, units=["mV"]*n_sig,
                    sig_name=["ch1otv1", "ch1otv2", "ch1otv3", "ch2ovt1", "ch2otv2", "ch2otv3"],
                    p_signal = signal, fmt=["16"]*n_sig, write_dir="directory")
        signals, fields = wfdb.rdsamp("directory\name", sampfrom=0, sampto="end", channels="all")
        
        fs = 4000
        t = np.arange(0, len(signals)/fs, 1.0/fs)
        plt.plot(t, signals)
        plt.xlim([0, 0.01])
        plt.title("Signals")
        print("Look at the plots")
        plt.show()
        
    else:
        print("Задайте файл")
