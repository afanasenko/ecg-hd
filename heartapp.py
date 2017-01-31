import sys
import os

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

        #self.sc = MyStaticMplCanvas(self.ui.plotArea)

        self.ui.fileList.itemClicked.connect(self.chooseFile)


    def listFiles(self, dirname):

        self.cur_dir = dirname
        self.ui.fileList.clear()

        for f in os.listdir(dirname):
            fullfile = os.path.join(dirname, f)
            if os.path.isfile(fullfile):
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
            print(full_file_name)
            self.ui.plotArea.compute_initial_figure()

    # Реакция на выбор каталога
    def changeCurrentDirectory(self):
        dirname = QFileDialog.getExistingDirectory()
        if dirname is not None:
            self.listFiles(dirname)


if __name__ == "__main__":

    qApp = QtGui.QApplication(sys.argv)
    wnd = ApplicationWindow()
    wnd.show()
    sys.exit(qApp.exec())

    # пример генерации кода из ui файла
    # ~/anaconda3/bin/pyuic4 ecgui.ui > ui_ecg_mainwindow.py
