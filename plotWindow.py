# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'plotWindow.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets
import numpy as np
from scipy import signal
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT  as  NavigationToolbar)
from matplotlib import cm


class Ui_plotWindow(object):
    def __init__(self,data):
        self.data=data

    def test(self):
        #y = np.arange(0, 100, 1)
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.plot(self.data)
        self.MplWidget.canvas.draw_idle()

    def setup_parameters(self):
        # self.data = self.data_raw
        self.time = int(self.lineEdit_time.text())
        # self.ch1 = int(self.lineEdit_ch1.text()) - 1
        # self.ch2 = int(self.lineEdit_ch2.text()) - 1
        # input channel
        ch = self.lineEdit_ch.text()
        self.channel = np.empty(0, dtype=int)
        ch = ch.split(',')
        for i in range(len(ch)):
            if '-' in ch[i]:
                tmp = ch[i].split('-')
                tmp = range(int(tmp[0]), int(tmp[1]))
                self.channel = np.append(self.channel, tmp)
            else:
                self.channel = np.append(self.channel, int(ch[i]))

        self.fs = int(self.lineEdit_fs.text())
        self.label = np.zeros(60)
        self.cmax = 100

        # self.section=self.section[:,3000 * self.time: 3000 * (self.time + 1)]

        if self.checkBox_normalize.isChecked():
            for i in self.channel:
                self.data[i] = self.data[i] / np.max(self.data[i])
            self.cmax = 1

        self.data = self.data[self.channel]
        # QtWidgets.QMessageBox.about(self, '', 'Complete')

    def BPfilter(self):
        fs = self.fs
        f1 = float(self.lineEdit_freq1.text())
        f2 = float(self.lineEdit_freq2.text())
        nmin = self.data.shape[1] / (fs * 60)
        for i in range(self.data.shape[0]):
            trace0 = self.data[i]
            for j in range(int(nmin)):
                trace = trace0[j * fs * 60:fs * 60 + j * fs * 60]
                trace = trace - trace.mean()
                trace = signal.detrend(trace)
                trace0[j * fs * 60:fs * 60 + j * fs * 60] = trace
            b, a = signal.butter(6, [2 * f1 / fs, 2 * f2 / fs], 'bandpass')
            filtedData = signal.filtfilt(b, a, trace0)
            # filtedData = filtedData * np.hanning(trace0.shape[0])
            self.data[i] = filtedData

    def selectionchange(self, index):
        self.plot_type = self.cb_PlotType.itemText(index)

    def update_graph(self):
        nt = self.fs * 60
        # t=int(self.lineEdit_time.text())
        # ch1=int(self.lineEdit_ch1.text())-1
        # ch2=int(self.lineEdit_ch2.text())-1
        # section = self.data[self.channel,3000*self.time:3000*(self.time+1)]
        self.section = self.data[:, nt * self.time: nt * (self.time + 1)]
        self.plot_core()

    def plot_next(self):
        nt = self.fs * 60
        if self.time == 59:
            QtWidgets.QMessageBox.about(self, 'Error', 'This is the last slice')
        else:
            self.time = self.time + 1
            self.section = self.data[:, nt * self.time:nt * (self.time + 1)]
            self.plot_core()

    def plot_back(self):
        nt = self.fs * 60
        if self.time == 0:
            QtWidgets.QMessageBox.about(self, 'Error', 'This is the first slice')
        else:
            self.time = self.time - 1
            # section = self.data[self.channel,3000*self.time:3000*(self.time+1)]
            self.section = self.data[:, nt * self.time: nt * (self.time + 1)]
            self.plot_core()

    def test(self):
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.plot(self.ui.pd_series.pd_series.dat.T[967])
        self.MplWidget.canvas.draw_idle()

    def plot_core(self):
        nt = self.fs * 60
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.clear()
        self.plot_type = 'Color Plot'
        if self.plot_type == 'Color Plot':
            self.MplWidget.canvas.axes.imshow(self.section, vmin=-self.cmax, vmax=self.cmax, aspect='auto',
                                              cmap=cm.seismic)
            self.MplWidget.canvas.axes.xaxis.set(ticks=np.linspace(0, nt, 13),
                                                 ticklabels=np.linspace(0, 60, 13, dtype=int))
            self.MplWidget.canvas.axes.set(ylabel='Channel', xlabel='Time(s)')
            self.MplWidget.canvas.draw_idle()
        elif self.plot_type == 'Wiggle Plot':
            for i in range(self.section.shape[0]):
                times = np.linspace(0, 60, nt)
                self.MplWidget.canvas.axes.plot(times, self.section[i] + i, c="k", linewidth=0.5)
            self.MplWidget.canvas.axes.set(ylabel='Channel', xlabel='Time(s)')
            self.MplWidget.canvas.draw_idle()
        elif self.plot_type == 'Check NCFs':
            times = np.linspace(0, 60, nt)
            tr1 = self.section[0]
            tr2 = self.section[1]
            self.MplWidget.canvas.axes.plot(times, tr1 + 1, c="k", linewidth=0.5)
            self.MplWidget.canvas.axes.plot(times, tr2, c="k", linewidth=0.5)
            self.MplWidget.canvas.draw_idle()
            times0 = np.linspace(-5, 5, 10 * self.fs + 1)
            self.MplWidget.canvas.axes.plot(times0, 1000 * self.xcorr[self.time], c="k", linewidth=0.5)
            self.MplWidget.canvas.draw_idle()

    def start_pick(self):
        self.label = np.zeros(60)

    def pick(self):
        self.label[self.time] = 1

    def reject(self):
        self.label[self.time] = 0

    def export(self):
        np.save('/Users/shenjunzhu/DAS/data/label', self.label)
        QtWidgets.QMessageBox.about(self, '', 'success')



    def setupUi(self, plotWindow):
        plotWindow.setObjectName("plotWindow")
        plotWindow.resize(1125, 700)
        self.centralwidget = QtWidgets.QWidget(plotWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayoutWidget_2 = QtWidgets.QWidget(self.centralwidget)
        self.horizontalLayoutWidget_2.setGeometry(QtCore.QRect(70, 70, 141, 32))
        self.horizontalLayoutWidget_2.setObjectName("horizontalLayoutWidget_2")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_2)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.pushButton_plot = QtWidgets.QPushButton(self.horizontalLayoutWidget_2)
        self.pushButton_plot.setObjectName("pushButton_plot")
        self.horizontalLayout_2.addWidget(self.pushButton_plot)
        self.MplWidget = MplWidget(self.centralwidget)
        self.MplWidget.setGeometry(QtCore.QRect(240, 60, 861, 581))
        self.MplWidget.setMinimumSize(QtCore.QSize(480, 320))
        self.MplWidget.setObjectName("MplWidget")
        self.lineEdit_time = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_time.setGeometry(QtCore.QRect(100, 180, 111, 21))
        self.lineEdit_time.setTabletTracking(False)
        self.lineEdit_time.setObjectName("lineEdit_time")
        self.lineEdit_ch = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_ch.setGeometry(QtCore.QRect(100, 210, 111, 21))
        self.lineEdit_ch.setObjectName("lineEdit_ch")
        self.lineEdit_freq1 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_freq1.setGeometry(QtCore.QRect(110, 290, 51, 21))
        self.lineEdit_freq1.setObjectName("lineEdit_freq1")
        self.lineEdit_freq2 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_freq2.setGeometry(QtCore.QRect(170, 290, 51, 21))
        self.lineEdit_freq2.setObjectName("lineEdit_freq2")
        self.pushButton_back = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_back.setGeometry(QtCore.QRect(230, 0, 141, 31))
        self.pushButton_back.setObjectName("pushButton_back")
        self.pushButton_next = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_next.setGeometry(QtCore.QRect(960, 0, 141, 31))
        self.pushButton_next.setObjectName("pushButton_next")
        self.pushButton_reject = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_reject.setGeometry(QtCore.QRect(230, 30, 141, 31))
        self.pushButton_reject.setObjectName("pushButton_reject")
        self.pushButton_pick = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_pick.setGeometry(QtCore.QRect(960, 30, 141, 31))
        self.pushButton_pick.setObjectName("pushButton_pick")
        self.pushButton_export = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_export.setGeometry(QtCore.QRect(60, 480, 141, 31))
        self.pushButton_export.setObjectName("pushButton_export")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(10, 180, 59, 16))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(10, 210, 59, 16))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(10, 290, 91, 16))
        self.label_3.setObjectName("label_3")
        self.checkBox_normalize = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_normalize.setGeometry(QtCore.QRect(70, 250, 111, 20))
        self.checkBox_normalize.setObjectName("checkBox_normalize")
        self.cb_PlotType = QtWidgets.QComboBox(self.centralwidget)
        self.cb_PlotType.setGeometry(QtCore.QRect(80, 110, 121, 32))
        self.cb_PlotType.setObjectName("cb_PlotType")
        self.cb_PlotType.addItem("")
        self.cb_PlotType.addItem("")
        self.pushButton_preprocess = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_preprocess.setGeometry(QtCore.QRect(80, 320, 113, 32))
        self.pushButton_preprocess.setObjectName("pushButton_preprocess")
        self.lineEdit_fs = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_fs.setGeometry(QtCore.QRect(100, 150, 51, 21))
        self.lineEdit_fs.setObjectName("lineEdit_fs")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(10, 150, 59, 16))
        self.label_4.setObjectName("label_4")
        self.pushButton_check = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_check.setGeometry(QtCore.QRect(80, 380, 112, 32))
        self.pushButton_check.setObjectName("pushButton_check")
        plotWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(plotWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1125, 22))
        self.menubar.setObjectName("menubar")
        plotWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(plotWindow)
        self.statusbar.setObjectName("statusbar")
        plotWindow.setStatusBar(self.statusbar)


        self.pushButton_plot.clicked.connect(self.update_graph)
        self.pushButton_check.clicked.connect(self.setup_parameters)
        self.pushButton_back.clicked.connect(self.plot_back)
        self.pushButton_next.clicked.connect(self.plot_next)
        self.pushButton_pick.clicked.connect(self.pick)
        self.pushButton_reject.clicked.connect(self.reject)
        self.pushButton_export.clicked.connect(self.export)
        self.pushButton_preprocess.clicked.connect(self.BPfilter)

        self.cb_PlotType.currentIndexChanged.connect(self.selectionchange)

        self.retranslateUi(plotWindow)
        QtCore.QMetaObject.connectSlotsByName(plotWindow)

    def retranslateUi(self, plotWindow):
        _translate = QtCore.QCoreApplication.translate
        plotWindow.setWindowTitle(_translate("plotWindow", "MainWindow"))
        self.pushButton_plot.setText(_translate("plotWindow", "Plot"))
        self.lineEdit_time.setText(_translate("plotWindow", "0"))
        self.lineEdit_ch.setText(_translate("plotWindow", "0-2120"))
        self.lineEdit_freq1.setText(_translate("plotWindow", "0.5"))
        self.lineEdit_freq2.setText(_translate("plotWindow", "24"))
        self.pushButton_back.setText(_translate("plotWindow", "Back"))
        self.pushButton_next.setText(_translate("plotWindow", "Next"))
        self.pushButton_reject.setText(_translate("plotWindow", "Reject"))
        self.pushButton_pick.setText(_translate("plotWindow", "Pick"))
        self.pushButton_export.setText(_translate("plotWindow", "Export"))
        self.label.setText(_translate("plotWindow", "Start_min"))
        self.label_2.setText(_translate("plotWindow", "Channel"))
        self.label_3.setText(_translate("plotWindow", "Bandpass (Hz)"))
        self.checkBox_normalize.setText(_translate("plotWindow", "Normalization"))
        self.cb_PlotType.setItemText(0, _translate("plotWindow", "Color Plot"))
        self.cb_PlotType.setItemText(1, _translate("plotWindow", "Wiggle Plot"))
        self.pushButton_preprocess.setText(_translate("plotWindow", "Preprocess"))
        self.lineEdit_fs.setText(_translate("plotWindow", "500"))
        self.label_4.setText(_translate("plotWindow", "fs"))
        self.pushButton_check.setText(_translate("plotWindow", "Confirm"))
from mplwidget import MplWidget

#     def setupUi(self, plotWindow):
#         plotWindow.setObjectName("plotWindow")
#         plotWindow.resize(1125, 700)
#         self.centralwidget = QtWidgets.QWidget(plotWindow)
#         self.centralwidget.setObjectName("centralwidget")
#         self.horizontalLayoutWidget_2 = QtWidgets.QWidget(self.centralwidget)
#         self.horizontalLayoutWidget_2.setGeometry(QtCore.QRect(70, 70, 141, 32))
#         self.horizontalLayoutWidget_2.setObjectName("horizontalLayoutWidget_2")
#         self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_2)
#         self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
#         self.horizontalLayout_2.setObjectName("horizontalLayout_2")
#         self.pushButton_plot = QtWidgets.QPushButton(self.horizontalLayoutWidget_2)
#         self.pushButton_plot.setObjectName("pushButton_plot")
#         self.horizontalLayout_2.addWidget(self.pushButton_plot)
#         self.MplWidget = MplWidget(self.centralwidget)
#         self.MplWidget.setGeometry(QtCore.QRect(240, 60, 861, 581))
#         self.MplWidget.setMinimumSize(QtCore.QSize(480, 320))
#         self.MplWidget.setObjectName("MplWidget")
#         self.lineEdit_time = QtWidgets.QLineEdit(self.centralwidget)
#         self.lineEdit_time.setGeometry(QtCore.QRect(70, 180, 111, 21))
#         self.lineEdit_time.setTabletTracking(False)
#         self.lineEdit_time.setObjectName("lineEdit_time")
#         self.lineEdit_ch = QtWidgets.QLineEdit(self.centralwidget)
#         self.lineEdit_ch.setGeometry(QtCore.QRect(70, 210, 111, 21))
#         self.lineEdit_ch.setObjectName("lineEdit_ch")
#         self.lineEdit_freq1 = QtWidgets.QLineEdit(self.centralwidget)
#         self.lineEdit_freq1.setGeometry(QtCore.QRect(70, 360, 51, 21))
#         self.lineEdit_freq1.setObjectName("lineEdit_freq1")
#         self.lineEdit_freq2 = QtWidgets.QLineEdit(self.centralwidget)
#         self.lineEdit_freq2.setGeometry(QtCore.QRect(130, 360, 51, 21))
#         self.lineEdit_freq2.setObjectName("lineEdit_freq2")
#         self.pushButton_back = QtWidgets.QPushButton(self.centralwidget)
#         self.pushButton_back.setGeometry(QtCore.QRect(230, 0, 141, 31))
#         self.pushButton_back.setObjectName("pushButton_back")
#         self.pushButton_next = QtWidgets.QPushButton(self.centralwidget)
#         self.pushButton_next.setGeometry(QtCore.QRect(960, 0, 141, 31))
#         self.pushButton_next.setObjectName("pushButton_next")
#         self.pushButton_reject = QtWidgets.QPushButton(self.centralwidget)
#         self.pushButton_reject.setGeometry(QtCore.QRect(230, 30, 141, 31))
#         self.pushButton_reject.setObjectName("pushButton_reject")
#         self.pushButton_pick = QtWidgets.QPushButton(self.centralwidget)
#         self.pushButton_pick.setGeometry(QtCore.QRect(960, 30, 141, 31))
#         self.pushButton_pick.setObjectName("pushButton_pick")
#         self.pushButton_export = QtWidgets.QPushButton(self.centralwidget)
#         self.pushButton_export.setGeometry(QtCore.QRect(60, 480, 141, 31))
#         self.pushButton_export.setObjectName("pushButton_export")
#         self.label = QtWidgets.QLabel(self.centralwidget)
#         self.label.setGeometry(QtCore.QRect(10, 180, 59, 16))
#         self.label.setObjectName("label")
#         self.label_2 = QtWidgets.QLabel(self.centralwidget)
#         self.label_2.setGeometry(QtCore.QRect(10, 210, 59, 16))
#         self.label_2.setObjectName("label_2")
#         self.label_3 = QtWidgets.QLabel(self.centralwidget)
#         self.label_3.setGeometry(QtCore.QRect(10, 360, 59, 16))
#         self.label_3.setObjectName("label_3")
#         self.checkBox_normalize = QtWidgets.QCheckBox(self.centralwidget)
#         self.checkBox_normalize.setGeometry(QtCore.QRect(70, 250, 111, 20))
#         self.checkBox_normalize.setObjectName("checkBox_normalize")
#         self.cb_PlotType = QtWidgets.QComboBox(self.centralwidget)
#         self.cb_PlotType.setGeometry(QtCore.QRect(70, 110, 121, 32))
#         self.cb_PlotType.setObjectName("cb_PlotType")
#         self.cb_PlotType.addItem("")
#         self.cb_PlotType.addItem("")
#         self.pushButton_Bandpass = QtWidgets.QPushButton(self.centralwidget)
#         self.pushButton_Bandpass.setGeometry(QtCore.QRect(70, 400, 113, 32))
#         self.pushButton_Bandpass.setObjectName("pushButton_Bandpass")
#         self.lineEdit_fs = QtWidgets.QLineEdit(self.centralwidget)
#         self.lineEdit_fs.setGeometry(QtCore.QRect(70, 150, 51, 21))
#         self.lineEdit_fs.setObjectName("lineEdit_fs")
#         self.label_4 = QtWidgets.QLabel(self.centralwidget)
#         self.label_4.setGeometry(QtCore.QRect(10, 150, 59, 16))
#         self.label_4.setObjectName("label_4")
#         self.pushButton_test = QtWidgets.QPushButton(self.centralwidget)
#         self.pushButton_test.setGeometry(QtCore.QRect(60, 440, 113, 32))
#         self.pushButton_test.setObjectName("pushButton_test")
#         plotWindow.setCentralWidget(self.centralwidget)
#         self.menubar = QtWidgets.QMenuBar(plotWindow)
#         self.menubar.setGeometry(QtCore.QRect(0, 0, 1125, 22))
#         self.menubar.setObjectName("menubar")
#         plotWindow.setMenuBar(self.menubar)
#         self.statusbar = QtWidgets.QStatusBar(plotWindow)
#         self.statusbar.setObjectName("statusbar")
#         plotWindow.setStatusBar(self.statusbar)
#
#         self.retranslateUi(plotWindow)
#         QtCore.QMetaObject.connectSlotsByName(plotWindow)
#
#         self.pushButton_test.clicked.connect(self.test)
#
#     def retranslateUi(self, plotWindow):
#         _translate = QtCore.QCoreApplication.translate
#         plotWindow.setWindowTitle(_translate("plotWindow", "MainWindow"))
#         self.pushButton_plot.setText(_translate("plotWindow", "Plot"))
#         self.lineEdit_time.setText(_translate("plotWindow", "0"))
#         self.lineEdit_ch.setText(_translate("plotWindow", "0-2120"))
#         self.lineEdit_freq1.setText(_translate("plotWindow", "0.5"))
#         self.lineEdit_freq2.setText(_translate("plotWindow", "24"))
#         self.pushButton_back.setText(_translate("plotWindow", "Back"))
#         self.pushButton_next.setText(_translate("plotWindow", "Next"))
#         self.pushButton_reject.setText(_translate("plotWindow", "Reject"))
#         self.pushButton_pick.setText(_translate("plotWindow", "Pick"))
#         self.pushButton_export.setText(_translate("plotWindow", "Export"))
#         self.label.setText(_translate("plotWindow", "Minute"))
#         self.label_2.setText(_translate("plotWindow", "Channel"))
#         self.label_3.setText(_translate("plotWindow", "Freq"))
#         self.checkBox_normalize.setText(_translate("plotWindow", "Normalization"))
#         self.cb_PlotType.setItemText(0, _translate("plotWindow", "Color Plot"))
#         self.cb_PlotType.setItemText(1, _translate("plotWindow", "Wiggle Plot"))
#         self.pushButton_Bandpass.setText(_translate("plotWindow", "BandPass"))
#         self.lineEdit_fs.setText(_translate("plotWindow", "500"))
#         self.label_4.setText(_translate("plotWindow", "Fs"))
#         self.pushButton_test.setText(_translate("plotWindow", "test"))
#
# from mplwidget import MplWidget