# coding: utf-8

import io
import os
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
    """
    Чтение записи zetlab, состоящей из тре файлов - .ana (данные),
    .anp (заголовок), .xml (тоже заголовок)
    :param f: полное или базовое (без расширения) имя файла
    :return: массив отсчетов (в физ. единицах), частота дискретизации, метка
    """

    if os.path.isfile(f):
        basename = ".".join(f.split(".")[:-1])
    else:
        basename = f

    if all([os.path.isfile(basename+x) for x in (".xml", ".ana", ".anp")]):
        pass
    else:
        print("Отсутствуют необходимые файлы")

    kv = load_xml_header(basename)

    if not kv:
        #anp может содержать кириллицу в cp1251
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


if __name__ == "__main__":

    fn = "/Users/arseniy/GUAP-ZAYC/rabbit20171013/s171012_141529/sig0003.ana"
    d, fs, n = anaread(fn)

    t = np.arange(0, len(d)/fs, 1.0/fs)

    plt.plot(t,d)
    plt.xlim([540, 545])
    plt.title(n)

    print("Look at the plots")
    plt.show()

