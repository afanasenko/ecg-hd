#coding: utf-8

import io
import os
import array
from matplotlib import pyplot as plt
import numpy as np
import xml.etree.ElementTree as ET


def load_xml_header(basename):
    kv = {}
    xmlname = basename+".xml"
    if os.path.isfile(xmlname):
        tree = ET.parse(xmlname)
        root = tree.getroot()
        desc = root.find("Signal")
        if desc is not None:
            kv["FRQ"] = int(desc.get("frequency",1))
            kv["ABSVOLT"] = float(desc.get("minlevel", 1))
            kv["CONVERT"] = desc.get("conversion", "mV")
            kv["NAME"] = desc.get("name", "unknown")

    return kv


def readbin(f):

    kv = load_xml_header(f)

    if not kv:
        with open(f + ".anp", "r") as fh:
            for line in fh:
                if line.strip():
                    p = [x.strip() for x in line.split()]
                    if len(p) == 2:
                        kv[p[0]] = p[1]

    data = array.array('f')

    with open(f+".ana", "rb") as fh:
        fh.seek(0, io.SEEK_END)
        filesz = fh.tell()
        num_samples = int(filesz/4)
        fh.seek(0, io.SEEK_SET)
        data.fromfile(fh, num_samples)

    fs = float(kv["FRQ"])
    scale = float(kv["ABSVOLT"])
    if kv["CONVERT"] == "mV":
        scale *= 1000

    return scale*np.array(data), fs, kv.get("name", "unknown")

if __name__ == "__main__":
    fn = "/Users/arseniy/Downloads/кроль_13_окт/s171012_141529/sig0003"
    #fn = "/Users/arseniy/GUAP-ZAYC/rabbit20170602/Кролик 02.06.2017 операция/sig0001"
    d, fs, n = readbin(fn)

    t = np.arange(0, len(d)/fs, 1.0/fs)

    plt.plot(t,d)
    plt.xlim([540, 545])
    plt.title(n)
    print("Look at the plots")
    plt.savefig("отведение 3 тромб через 2 мин.png")
