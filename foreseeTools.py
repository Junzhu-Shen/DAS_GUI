# Last edit: Oct 28, 2020

import numpy as np
from scipy import signal
import scipy.fftpack
import datetime
import matplotlib.pyplot as plt

def downsample(data,rate=10,axis=1):
    length=data.shape[axis]
    index=np.arange(0,length,int(rate))
    if axis==0:
        data_downsampled=data[index]
    if axis==1:
        data_downsampled=data[:,index]
    return data_downsampled

def bpfilter(trace,f1,f2,fs):
    b,a=signal.butter(4,[2*f1/fs,2*f2/fs],'bandpass')
    trace=signal.filtfilt(b,a,trace)
    return trace

def freqSpectr(trace,fs):
    N=trace.size
    xf=np.linspace(0.0,fs/2.0,N//2)
    yf=scipy.fftpack.fft(trace)
    plt.plot(xf,abs(yf[0:N//2]))
    plt.show()

def read_raw(stime, etime, fs, datadir, skip = 1, channelFile='./DAS_location.txt'):
    from tdms_reader import TdmsReader
    location=np.loadtxt(channelFile)
    ch_info=location[:,1].astype(int)
    n_file=int((etime-stime).total_seconds() / 60.0)
    data=np.empty([2120//skip,30000*n_file])
    for i in range(n_file):
        timestamp=str(stime+datetime.timedelta(minutes=i))
        timestamp=timestamp.replace('-','').replace(' ','_').replace(':','')[0:-3]
        inputFile=datadir+'PSUDAS_UTC_'+timestamp+'.tdms'
        tdms = TdmsReader(inputFile)
        props = tdms.get_properties()
        nCh = tdms.fileinfo['n_channels']
        # nSamples = tdms.channel_length
        if nCh==2432:
            data0 = tdms.get_data(0,nCh,0,tdms.channel_length-1)
        else:
            data0 = tdms.get_data(64,nCh,0,tdms.channel_length-1)
        # kill traces
        data0=np.transpose(data0)[ch_info]
        # spatial downsample
        if skip != 1:
            data0=data0[np.arange(0,2120,int(skip))]

        data[:,30000*i:30000*(i+1)]=data0
    data=downsample(data,500/fs)
    return data

### PSD ###
def computePSD(trace,param,fs,smooth=True):
    import math
    delta=1/fs
    winLen=param['winLen']
    nSeg=int(param['trLen']/winLen)
    nShift=param['winShift']*fs

    nSamp=2**int(math.log(winLen*fs,2))
    taperWindow=np.hanning(nSamp)
    startIndex=0
    m11 = np.zeros((nSamp//2)+1,dtype=np.complex)
    for n in range(nSeg):
        endIndex=startIndex+nSamp
        channelSegment=trace[startIndex:endIndex]
        channelSegment=channelSegment-np.mean(channelSegment)
        channelSegment=channelSegment * taperWindow
        FFT = np.zeros((nSamp//2)+1, dtype=np.complex)
        FFT = np.fft.rfft(channelSegment)
        m11 += FFT[0:(nSamp//2)+1] * np.conjugate(FFT[0:(nSamp//2)+1])
        startIndex += nShift
    # Convert FFT to PSD
    norm  = 2.0 * delta /  float(nSamp)
    m11 /= float(nSeg)
    power = norm * np.abs(m11[0:(nSamp//2)+1])

    smoothX    = []
    smoothPSD  = []
    # Smooth
    if smooth:
        xType = 'frequency'
        maxT=10
        octaveWindowWidth=1.0/2.0
        octaveWindowShift=1.0/8.0
        maxPeriod=maxT * pow(2,octaveWindowWidth/2.0)    # maximum period needed to compute value at maxT period point
        #minFrequency   = 1.0/float(maxPeriod)
        minFrequency   = 0.05

        frequency = np.array( np.arange(0,(nSamp/2)+1)/float(nSamp * delta))

        smoothX,smoothPSD = smoothNyquist(xType,frequency,power,fs,
                                              octaveWindowWidth, octaveWindowShift, minFrequency)
        # smoothX,smoothPSD       = SFL.smoothF(frequency,power,fs, octaveWindowWidth, octaveWindowShift, minFrequency,25)

        # smoothPSD = 10.0* np.log10(smoothPSD)
    else:
        # smoothX = np.linspace(0,fs/2,nSamp//2)
        smoothX = np.fft.fftfreq(nSamp,d=delta)[slice(1,nSamp//2)]
        smoothPSD = power[:-2]
        # smoothPSD = 10*np.log10(power)
    return smoothX, smoothPSD

def getBin(X,Y,Xc,octaveHalfWindow):
    import math
    thisBin = []
    shift = math.pow(2.0,octaveHalfWindow)
    #
    # the bin is octaveHalfWindow around Xc
    #
    X1 = Xc / shift
    X2 = Xc * shift
    Xs = min(X1,X2)
    Xe = max(X1,X2)

    #
    # gather the values that fall within the range >=Xs and <= Xe
    #
    for i in range(0,len(X)):
        if X[i] >= Xs and X[i] <= Xe:
            thisBin.append(float(Y[i]))

    return thisBin

def smoothNyquist(type,Xi,Yi,samplingRate,octaveWindowWidth,octaveWindowShift,xLimit):
    import math
    X  = []
    Y  = []
    #
    # shortest period (highest frequency) center of the window
    # starts  at the Nyquist
    #
    windowWidth = octaveWindowWidth
    halfWindow  = float(windowWidth / 2.0)
    windowShift = octaveWindowShift
    shift       = math.pow(2.0,windowShift)        # shift of each Fc

    #
    # the first center X at the Nyquist
    #
    if type == "frequency":
        Xc = float(samplingRate) / float(2.0) # Nyquist frequency
    else:
        Xc = float(2.0) /float(samplingRate)  # Nyquist period

    while (True):
        #
        # do not go below the minimum frequency
        # do not go above the maximum period
        #
        if (type == "frequency" and Xc < xLimit) or (type == "period" and Xc > xLimit):
            break

        thisBin = getBin(Xi,Yi,Xc,halfWindow)

        #
        # bin should not be empty
        #
        if (len(thisBin) > 0):
            Y.append(np.mean(thisBin))
            X.append(Xc)
        else:
            Y.append(float('NAN'))
            X.append(Xc)

        #
        # move the center frequency to the right by half of the windowWidth
        # move the center period to the left by half of the windowWidth
        #
        if type == "frequency":
            Xc /= shift
        else:
            Xc *= shift
    #
    # sort on X and return
    #
    X,Y = (list(t) for t in zip(*sorted(zip(X,Y))))
    return (X,Y)

### ### ###

def preprocess(data, f1, f2, fs, w=125, order=4, bandpass=True, onebit=False, whitening=False):
    # f1-f2: bandpass frequency range
    # w: set w=0 for 1-bit
    length = data.shape[1]
    nch = data.shape[0]
    data_new = np.zeros_like(data)
    b, a = signal.butter(order, [2 * f1 / fs, 2 * f2 / fs], 'bandpass')
    for i in range(nch):
        # demean/detrend
        trace = data[i]
        trace = trace - trace.mean()
        trace = signal.detrend(trace)
        if bandpass:
            trace = signal.filtfilt(b, a, trace)
        # running-absolute-mean normalization
        if onebit:
            trace=trace/np.abs(trace)

        # whitening
        if whitening:
            x = whiten(trace, length, 1 / fs, f1, f2, plot=False)
            trace = np.real(scipy.fftpack.ifft(x))

        data_new[i] = trace
    return data_new

def pre_1st(data, f1=0.5, f2=24, f1_0=5, f2_0=25, fs=50, w=125, order=4, bandpass=True, RAM=True,
            bp_before_RAM=False, median_filter=False, whitening=True):
    # f1-f2: bandpass frequency range
    # w: set w=0 for 1-bit
    length = data.shape[1]
    nch = data.shape[0]
    data_new = np.zeros_like(data)
    for i in range(nch):
        # demean/detrend
        trace = data[i]
        trace = trace - trace.mean()
        trace = signal.detrend(trace)
        if bandpass:
            b, a = signal.butter(order, [2 * f1 / fs, 2 * f2 / fs], 'bandpass')
            trace = signal.filtfilt(b, a, trace)
        # running-absolute-mean normalization
        if RAM:
            trace=trace/np.abs(trace)
            # weight = np.zeros_like(trace)
            # if bp_before_RAM:
            #     b, a = signal.butter(order, [2 * f1_0 / fs, 2 * f2_0 / fs], 'bandpass')
            #     filtedData = signal.filtfilt(b, a, trace)
            #     trace0 = np.zeros(length + 2 * w)
            #     trace0[w:w + length] = filtedData
            # else:
            #     trace0 = np.zeros(length + 2 * w)
            #     trace0[w:w + length] = trace
            # for j in range(length):
            #     wintrace = trace0[j:j + 2 * w + 1]
            #     weight[j] = np.sum(np.abs(wintrace)) / (2 * w + 1)
            # trace = trace / weight

        ####

        if median_filter:
            trace = scipy.signal.medfilt(trace, 51)

        # whitening
        if whitening:
            x = whiten(trace, length, 1 / fs, f1, f2, plot=False)
            trace = np.real(scipy.fftpack.ifft(x))

        data_new[i] = trace
    return data_new


def whiten(data, Nfft, delta, freqmin, freqmax, plot=False):
    """This function takes 1-dimensional *data* timeseries array,
    goes to frequency domain using fft, whitens the amplitude of the spectrum
    in frequency domain between *freqmin* and *freqmax*
    and returns the whitened fft.

    :type data: :class:`numpy.ndarray`
    :param data: Contains the 1D time series to whiten
    :type Nfft: int
    :param Nfft: The number of points to compute the FFT
    :type delta: float
    :param delta: The sampling frequency of the `data`
    :type freqmin: float
    :param freqmin: The lower frequency bound
    :type freqmax: float
    :param freqmax: The upper frequency bound
    :type plot: bool
    :param plot: Whether to show a raw plot of the action (default: False)

    :rtype: :class:`numpy.ndarray`
    :returns: The FFT of the input trace, whitened between the frequency bounds
"""

    if plot:
        plt.subplot(411)
        plt.plot(np.arange(len(data)) * delta, data)
        plt.xlim(0, len(data) * delta)
        plt.title('Input trace')

    Napod = 100
    Nfft = int(Nfft)
    freqVec = scipy.fftpack.fftfreq(Nfft, d=delta)[:Nfft // 2]

    J = np.where((freqVec >= freqmin) & (freqVec <= freqmax))[0]
    low = J[0] - Napod
    if low <= 0:
        low = 1

    porte1 = J[0]
    porte2 = J[-1]
    high = J[-1] + Napod
    if high > Nfft / 2:
        high = int(Nfft // 2)

    FFTRawSign = scipy.fftpack.fft(data, Nfft)

    if plot:
        plt.subplot(412)
        axis = np.arange(len(FFTRawSign))
        plt.plot(axis[1:], np.abs(FFTRawSign[1:]))
        plt.xlim(0, max(axis))
        plt.title('FFTRawSign')

    # Left tapering:
    FFTRawSign[0:low] *= 0
    FFTRawSign[low:porte1] = np.cos(
        np.linspace(np.pi / 2., np.pi, porte1 - low)) ** 2 * np.exp(
        1j * np.angle(FFTRawSign[low:porte1]))
    # Pass band:
    FFTRawSign[porte1:porte2] = np.exp(1j * np.angle(FFTRawSign[porte1:porte2]))
    # Right tapering:
    FFTRawSign[porte2:high] = np.cos(
        np.linspace(0., np.pi / 2., high - porte2)) ** 2 * np.exp(
        1j * np.angle(FFTRawSign[porte2:high]))
    FFTRawSign[high:Nfft + 1] *= 0

    # Hermitian symmetry (because the input is real)
    FFTRawSign[-(Nfft // 2) + 1:] = FFTRawSign[1:(Nfft // 2)].conjugate()[::-1]

    if plot:
        plt.subplot(413)
        axis = np.arange(len(FFTRawSign))
        plt.axvline(low, c='g')
        plt.axvline(porte1, c='g')
        plt.axvline(porte2, c='r')
        plt.axvline(high, c='r')

        plt.axvline(Nfft - high, c='r')
        plt.axvline(Nfft - porte2, c='r')
        plt.axvline(Nfft - porte1, c='g')
        plt.axvline(Nfft - low, c='g')

        plt.plot(axis, np.abs(FFTRawSign))
        plt.xlim(0, max(axis))

        wdata = np.real(scipy.fftpack.ifft(FFTRawSign, Nfft))
        plt.subplot(414)
        plt.plot(np.arange(len(wdata)) * delta, wdata)
        plt.xlim(0, len(wdata) * delta)
        plt.show()
    return FFTRawSign


def txcorr(x1, x2, lag):
    length = x1.size
    corr = np.correlate(x1, x2, 'full')
    xcorr = corr[length - 1 - lag:length + lag]
    return xcorr

# plot Function
def beamformingPlot(out):
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize

    from obspy import Trace,Stream
    from obspy.core.util import AttribDict
    from obspy import UTCDateTime
    from obspy.imaging.cm import obspy_sequential
    from obspy.signal.invsim import corn_freq_2_paz
    from obspy.signal.array_analysis import array_processing

    cmap = obspy_sequential

    # make output human readable, adjust backazimuth to values between 0 and 360
    t, rel_power, abs_power, baz, slow = out.T
    baz[baz < 0.0] += 360

    # choose number of fractions in plot (desirably 360 degree/N is an integer!)
    N = 36
    N2 = 30
    abins = np.arange(N + 1) * 360. / N
    sbins = np.linspace(0, 3, N2 + 1)

    # sum rel power in bins given by abins and sbins
    hist, baz_edges, sl_edges = \
        np.histogram2d(baz, slow, bins=[abins, sbins], weights=rel_power)

    # transform to radian
    baz_edges = np.radians(baz_edges)

    # add polar and colorbar axes
    fig = plt.figure(figsize=(8, 8))
    cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
    ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("N")

    dh = abs(sl_edges[1] - sl_edges[0])
    dw = abs(baz_edges[1] - baz_edges[0])

    # circle through backazimuth
    for i, row in enumerate(hist):
        bars = ax.bar((i * dw) * np.ones(N2),
                      height=dh * np.ones(N2),
                      width=dw, bottom=dh * np.arange(N2),
                      color=cmap(row / hist.max()))

    ax.set_xticks(np.linspace(0, 2 * np.pi, 4, endpoint=False))
    ax.set_xticklabels(['N', 'E', 'S', 'W'],fontsize=16)

    # set slowness limits
    ax.set_ylim(0, 3)
    [i.set_color('grey') for i in ax.get_yticklabels()]

    ColorbarBase(cax, cmap=cmap,
                 norm=Normalize(vmin=hist.min(), vmax=hist.max()))
