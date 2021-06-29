from PyQt5.uic import loadUi
from PyQt5 import QtCore, QtGui, QtWidgets
from plotWindow import Ui_plotWindow
from tdms_reader import TdmsReader
from foreseeTools import downsample

import datetime

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT  as  NavigationToolbar)
from matplotlib import cm

import numpy as np
import scipy.io

class Mymain(QtWidgets.QMainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        loadUi("LoadingWindow.ui", self)
        # self.ui=Ui_MainWindow()
        # self.ui.setupUi(self)

        self.setWindowTitle("DAS Data Visualization")

        # self.pushButton_load.clicked.connect(self.getCSV)
        # self.pushButton_plot.clicked.connect(self.update_graph)
        # self.pushButton_setup_parameters.clicked.connect(self.setup_parameters)
        # self.pushButton_back.clicked.connect(self.plot_back)
        # self.pushButton_next.clicked.connect(self.plot_next)
        # self.pushButton_pick.clicked.connect(self.pick)
        # self.pushButton_reject.clicked.connect(self.reject)
        # self.pushButton_export.clicked.connect(self.export)
        # self.pushButton_Bandpass.clicked.connect(self.BPfilter)
        # self.pushButton_loadtdms.clicked.connect(self.open2ndwindow)
        # self.pushButton_test.clicked.connect(self.test)

        self.pushButton_2ndwindow.clicked.connect(self.open2ndwindow)
        self.pushButton_start.clicked.connect(self.start)
        self.pushButton_getmat.clicked.connect(self.get_mat)

        self.data_raw=np.array
        self.label=np.array
        self.plot_type='Color Plot'

    def open2ndwindow(self):
        self.window = QtWidgets.QMainWindow()
        self.ui2=Ui_plotWindow(self.data)
        self.ui2.setupUi(self.window)
        self.window.addToolBar(NavigationToolbar(self.ui2.MplWidget.canvas, self))
        self.window.show()

    def start(self):
        datadir = '/Users/shenjunzhu/DAS/data/0830/'
        t1 = self.lineEdit_stime.text()
        t2 = self.lineEdit_etime.text()
        stime = datetime.datetime(int(t1.split(',')[0]),
                                  int(t1.split(',')[1]),
                                  int(t1.split(',')[2]),
                                  int(t1.split(',')[3]),
                                  int(t1.split(',')[4]),
                                  int(t1.split(',')[5]),
                                  int(t1.split(',')[6]) * 1000)
        etime = datetime.datetime(int(t2.split(',')[0]),
                                  int(t2.split(',')[1]),
                                  int(t2.split(',')[2]),
                                  int(t2.split(',')[3]),
                                  int(t2.split(',')[4]),
                                  int(t2.split(',')[5]),
                                  int(t2.split(',')[6]) * 1000)
        fs = int(self.lineEdit_fs.text())
        # for i in range(int((etime - stime).total_seconds() / 60.0)):
        #     timestamp = str(stime + datetime.timedelta(minutes=i))
        #     timestamp = timestamp.replace('-', '').replace(' ', '_').replace(':', '')[0:-3]
        #     file_tdms = datadir + 'PSUDAS_UTC_' + timestamp + '.tdms'
        #     pd_series = tdms2pd(file_tdms, traces=[0, 2432])
        #     self.data=np.arange(1,100)

        location = np.loadtxt('DAS_location.txt')
        ch_info = location[:, 1].astype(int)
        n_file = int((etime - stime).total_seconds() / 60.0)
        data = np.empty([2120, 30000 * n_file])
        for i in range(n_file):
            timestamp = str(stime + datetime.timedelta(minutes=i))
            timestamp = timestamp.replace('-', '').replace(' ', '_').replace(':', '')[0:-3]
            inputFile = datadir + 'PSUDAS_UTC_' + timestamp + '.tdms'
            tdms = TdmsReader(inputFile)
            props = tdms.get_properties()
            nCh = tdms.fileinfo['n_channels']
            # nSamples = tdms.channel_length
            if nCh == 2432:
                data0 = tdms.get_data(0, nCh, 0, tdms.channel_length - 1)
            else:
                data0 = tdms.get_data(64, nCh, 0, tdms.channel_length - 1)
            data[:, 30000 * i:30000 * (i + 1)] = np.transpose(data0)[ch_info]
        self.data = downsample(data, 500 / fs)

    def get_mat(self):
        filePath, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Select data", "/Users/shenjunzhu/DAS/data/0515", '*.mat')
        inputf = scipy.io.loadmat(str(filePath))
        self.data = np.array(inputf['data'])


if __name__=='__main__':
    app = QtWidgets.QApplication([])
    window = Mymain()
    window.show()
    app.exec_()
