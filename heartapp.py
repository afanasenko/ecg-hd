import sys
import os
import re
import wfdb
import json

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QFileDialog
from ui_mainwindow import Ui_MainWindow

"""
"""
class ApplicationWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.directorySelect.clicked.connect(self.changeCurrentDirectory)
        self.listFiles('.')

        self.loaded_signal = None
        self.loaded_signal_header = None

        self.ui.fileList.itemClicked.connect(self.chooseFile)
        self.ui.channelsDropdown.currentIndexChanged.connect(self.plotSignal)


    def listFiles(self, dirname):

        self.cur_dir = dirname
        self.ui.fileList.clear()

        for f in os.listdir(dirname):
            fullfile = os.path.join(dirname, f)
            # Не выводим каталоги
            if os.path.isfile(fullfile):
                # Выводим только .dat файлы
                file_re = re.compile(r'.*\.dat$')
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
        self.ui.samplingFreq.setText(str(fs))
        self.ui.signalUnits.setText(self.loaded_signal_header["units"][chan])

        samples = int(self.config_data["display"]["defaultTimeSpan"] * fs)
        # Если samples больше, чем размер файла - ошибки не будет
        self.ui.plotArea.plot_signal(self.loaded_signal[0:samples, chan], fs)

    def loadConfig(self, filename):
        self.config_data = {}
        with open(filename, "r") as fp:
            try:
                self.config_data = json.load(fp)
            except:
                print("Unable to load config")


if __name__ == "__main__":

    qApp = QtGui.QApplication(sys.argv)
    wnd = ApplicationWindow()
    wnd.loadConfig("config.json")

    wnd.show()
    sys.exit(qApp.exec())

    # пример генерации кода из ui файла
    # ~/anaconda3/bin/pyuic4 ecgui.ui > ui_ecg_mainwindow.py
