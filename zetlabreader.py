# coding: utf-8

import io
import os
import sys
import array
from matplotlib import pyplot as plt
import numpy as np
import xml.etree.ElementTree as ET


def load_xml_header(basename):
    kv = {}
    xmlname = basename + ".xml"
    if os.path.isfile(xmlname):
        tree = ET.parse(xmlname)
        root = tree.getroot()
        desc = root.find("Signal")
        if desc is not None:
            kv["FRQ"] = int(desc.get("frequency", 1))
            kv["ABSVOLT"] = float(desc.get("minlevel", 1))
            kv["CONVERT"] = desc.get("conversion", "mV")
            kv["NAME"] = desc.get("name", "unknown")

    return kv


def anaread(f):

    basename = ".".join(f.split(".")[:-1])

    kv = load_xml_header(basename)

    if not kv:
        with open(basename + ".anp", "r") as fh:
            for line in fh:
                if line.strip():
                    p = [x.strip() for x in line.split()]
                    if len(p) == 2:
                        kv[p[0]] = p[1]

    data = array.array('f')

    with open(basename+ ".ana", "rb") as fh:
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


def readint(f):
    b = f.read(4)
    i = int.from_bytes(b, byteorder="little", signed=True)
    return i


def readshort(f):
    b = f.read(2)
    i = int.from_bytes(b, byteorder="little", signed=True)
    return i


def holterread(filename):
    """
    read ecg data from .dat file in Holter format
    :param filename:
    :return:
    """
    with open(filename, "rb") as fh:
        marker = fh.read(8)
        crc = fh.read(2)
        varblock_len = readint(fh)
        ecg_samples = readint(fh)
        varblock_offset = readint(fh)
        ecg_offset = readint(fh)
        version = readshort(fh)


        print(varblock_len)
        print(ecg_samples)
        print(ecg_offset)

        ecg_offset=3940
        ecg_samples = 12000

        #data = array.array("i")  # signed int ???
        data = array.array("h")  # signed short
        #data = array.array("f")  # float
        fh.seek(ecg_offset, io.SEEK_SET)
        data.fromfile(fh, ecg_samples)

        print("...")
        raw = np.reshape(np.array(data), (-1,6))
        print(raw.shape)
        return raw[:, 2], 100, "lll"


if __name__ == "__main__":

    d, fs, n = holterread('/Users/arseniy/GUAP-ZAYC/holter/Чумичева.dat')

    #fn = "/Users/arseniy/Downloads/кроль_13_окт/s171012_141529/sig0003"
    #d, fs, n = anaread(fn)

    t = np.arange(0, len(d)/fs, 1.0/fs)


    plt.plot(t, d)
    #plt.xlim([540, 545])
    plt.title(n)

    print("Look at the plots")
    plt.show()
    #plt.savefig("отведение 3 тромб через 2 мин.png")
