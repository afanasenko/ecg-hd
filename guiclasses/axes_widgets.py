#!/usr/bin/env python
# coding: utf-8

from PyQt5 import QtWidgets

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

"""
    Обертка над одинарным графиком из matplotlib,
    унаследованная от QWidget
"""

class axesWidget(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        #
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(
            self,
            QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding
        )
        FigureCanvas.updateGeometry(self)

        self.compute_initial_figure()

    def compute_initial_figure(self):
        self.draw()

    def plot_signal(self, x, fs):
        t = np.array(range(0, x.shape[0]), float) / fs
        self.axes.plot(t, x)
        self.draw()
        pass

    def plot_signal_with_markup(self, x, marks, fs):
        t = np.array(range(0, x.shape[0]), float) / fs

        self.axes.plot(t, x, "b")
        self.axes.hold(True)
        self.axes.plot(marks/fs, [x[i] for i in marks], "r+", markersize=6)
        self.axes.hold(False)
        self.draw()
        pass