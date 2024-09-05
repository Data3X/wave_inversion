import obspy
import obspy.signal.filter
from obspy import read
import matplotlib.pyplot as plt
from math import *
import numpy as np
from scipy.fftpack import fft,ifft


filename = 'test1'
singlechannel = read('D:\\research\\dsm\\DSMsynTI-mpi-master\\results\\'+filename+'.Zs')

delta1 = singlechannel[0].stats.delta
data1 = singlechannel[0].data[0:8192]
T1 = [i*delta1 for i in range(len(data1))]

def read_sac(path):
    channel = read(path)
    data = channel[0]
    data.plot()
    plt.show()

def fourier_trans(T,sr):
    complex_array = np.fft.fft(sr)
    y_ = np.fft.ifft(complex_array).real
    Freqency = np.fft.fftfreq(y_.size,T[1]-T[0])
    Amplitude = np.abs(complex_array)
    Phase = np.angle(complex_array)
    return Freqency,Amplitude,Phase

def inverse_fourier_trans(Freqency,Amplitude):
    return 
    
def draw_Freqency_Amplitude(Freqency,Amplitude):
    plt.xscale("log")
    plt.yscale('log')
    plt.xlabel('Freqency/Hz')
    plt.ylabel('Amplitudelitude')
    plt.plot(Freqency[Freqency>0],Amplitude[Freqency>0])
    plt.show()

# plot seismogram, 'start_time' is a string, delta/s, x_unit:string, 
def draw_seismogram(signal,delta,x_unit='s',y_unit='m/s',title='Seismogram',start_time=None):
    time = [i*delta for i in range(len(signal))]
    plt.xlabel('Time'+'('+x_unit+')')
    plt.ylabel('Velocity'+'('+y_unit+')')
    plt.plot(time,signal)
    plt.show()

def filter(Freqency,Amplitude,center_period,half_period_width,attenuation):
    low_freq = 1/(center_period+half_period_width)
    high_freq = 1/(center_period-half_period_width)
    center_freq = 1/center_period
    for count in range(len(Freqency)):
        if Freqency[count] > low_freq and Freqency[count] < high_freq:
            Amplitude[count] = Amplitude[count]*2.718**(-attenuation*((Freqency[count]-center_freq)/center_freq)**2)
        elif Freqency[count] < -low_freq and Freqency[count] > -high_freq:
            Amplitude[count] = Amplitude[count]*2.718**(-attenuation*((Freqency[count]+center_freq)/center_freq)**2)
        else:
            Amplitude[count] = 0
    return Freqency,Amplitude


signal = read('seis_data\\XB.ELYSE.00.HHR.D.2022.124.230000.SAC')
data2 = signal[0].data
delta2 = signal[0].stats.delta
T2 = [i*delta2 for i in range(len(data2))]
def get_filtered_waveforms(T,data,delta):
    Freqency,Amplitude,Phase = fourier_trans(T,data)
    nperiod = 60
    arrival_time = []
    max_amplitude = []
    period_step = 0.5
    start_period = 10
    end_period = 60
    time_halfwidth = 5
    aut = 50.3
    nstep = round((end_period-start_period)/period_step)
    for point in range(nstep):
        freqency, amplitude = filter(Freqency,Amplitude/2,start_period+point*period_step,time_halfwidth,aut)
        freq_signal = amplitude*np.cos(Phase) + 1j*amplitude*np.sin(Phase)
        filtered_signal = np.real(np.fft.ifft(freq_signal))
        envelope = obspy.signal.filter.envelope(filtered_signal)
        enhance = 2*period_step
        # Normalization
        k = max(envelope)
        envelope1 = [enhance*i/k+start_period+point*period_step for i in envelope]
        max_amp = np.max(envelope)
        plt.plot([i*delta for i in range(len(envelope))],envelope1)
        max_amplitude.append(max_amp)
        arrival = np.argmax(envelope)*delta
        arrival_time.append(arrival)

    plt.xlabel('Time(s)')
    plt.ylabel('Period(s)')
    plt.title('Filterd Waveform')
    plt.show()

get_filtered_waveforms(T1,data1,delta1)


