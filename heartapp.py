#!/usr/bin/env python
# coding: utf-8

import sys
import os
import re
import wfdb
import json
import numpy as np

from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QApplication, QFileDialog, QMainWindow
from ui_mainwindow import Ui_MainWindow

from sigsegment import extract_short_peaks

"""
"""
class ApplicationWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.directorySelect.clicked.connect(self.changeCurrentDirectory)
        self.listFiles('.')

        self.loaded_signal = None
        self.loaded_signal_header = None

        self.ui.fileList.itemClicked.connect(self.chooseFile)
        self.ui.channelsDropdown.currentIndexChanged.connect(self.plotSignal)

        self.ui.refreshSignal.clicked.connect(self.plotSignal)

    def listFiles(self, dirname):

        self.cur_dir = dirname
        self.ui.fileList.clear()

        for f in os.listdir(dirname):
            fullfile = os.path.join(dirname, f)
            # Не выводим каталоги
            if os.path.isfile(fullfile):
                # Выводим только .hea файлы
                file_re = re.compile(r'.*\.hea$')
                if file_re.match(f) is not None:
                    self.ui.fileList.addItem(f)

        # Отображение имени текущего каталога на ярлычке
        max_len = 50
        if len(dirname) > max_len:
            dirname = '.../' + dirname.split(os.sep)[-1]
        self.ui.directoryLabel.setText(dirname)

    # Реакция на выбор файла
    def chooseFile(self, lwi):
        if lwi is not None:
            # Полное имя выбранного файла
            full_file_name = os.path.join(self.cur_dir, lwi.text())

            # Отрезаем расширение для чтения с помощью wfdb
            recordname = '.'.join(full_file_name.split('.')[:-1])

            # Чтение сигнала
            sig, fields = wfdb.rdsamp(recordname)
            # Не очень понятна обработка ошибок чтения

            self.loaded_signal = sig
            self.loaded_signal_header = fields

            self.updateChannelCombo(sig)
            self.plotSignal()

    # Реакция на выбор каталога
    def changeCurrentDirectory(self):
        dirname = QFileDialog.getExistingDirectory()
        if dirname and dirname is not None:
            self.listFiles(dirname)

    def updateChannelCombo(self, sig):
        numch = sig.shape[1]
        if self.ui.channelsDropdown.count() != numch:
            self.ui.channelsDropdown.clear()
            self.ui.channelsDropdown.addItems([str(x+1) for x in range(numch)])

    def plotSignal(self):
        chan = self.ui.channelsDropdown.currentIndex()
        fs = self.loaded_signal_header.get("fs", 360)
        self.ui.samplingFreq.setText(str(fs) + " Hz")

        dur = np.ceil(len(self.loaded_signal[:, chan]) / fs)
        self.ui.signalDuration.setText(str(dur) + " s")
        self.ui.signalUnits.setText(self.loaded_signal_header["units"][chan])

        samples = int(self.config_data["display"]["defaultTimeSpan"] * fs)
        # Если samples больше, чем размер файла - ошибки не будет

        if self.ui.autoProcessing.isChecked():
            # Считывание параметров привязки
            unbias_wnd = 0
            if self.ui.zeroBaseline.isChecked():
                try:
                    unbias_wnd = int(self.ui.biasWindowLen.text())
                except:
                    print("Unable to convert value")

            try:
                pk_interval = int(self.ui.minRinterval.text())
            except:
                pk_interval = self.config_data["rDetect"]["peakIntervalMs"]
                print("Unable to convert value")

            pk_len = self.config_data["rDetect"]["peakLengthMs"]

            pks, bks, hfs, lfs = extract_short_peaks(
                self.loaded_signal[0:samples, chan],
                fs,
                unbias_wnd,
                pk_len,
                pk_interval
            )

            if unbias_wnd:
                self.ui.plotArea.plot_signal_with_markup(
                    lfs, pks, fs
                )
            else:
                self.ui.plotArea.plot_signal_with_markup(
                    self.loaded_signal[0:samples, chan], pks, fs
                )
        else:
            self.ui.plotArea.plot_signal(self.loaded_signal[0:samples, chan], fs)

    def loadConfig(self, filename):
        self.config_data = {}
        with open(filename, "r") as fp:
            try:
                self.config_data = json.load(fp)
            except:
                print("Unable to load config")

        if self.config_data:
            self.ui.biasWindowLen.setText(
                str(self.config_data["preProcessing"]["biasWindowMs"])
            )
            self.ui.minRinterval.setText(
                str(self.config_data["rDetect"]["peakIntervalMs"])
            )

if __name__ == "__main__":


    app = QApplication(sys.argv)

    wnd = ApplicationWindow()
    wnd.loadConfig("config.json")

    wnd.show()
    sys.exit(app.exec_())

    # пример генерации кода из ui файла
    # ~/anaconda3/bin/pyuic4 ecgui.ui > ui_ecg_mainwindow.py
