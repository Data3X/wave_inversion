import obspy
import obspy.signal.filter
from obspy import read
import matplotlib.pyplot as plt
from math import *
import numpy as np
from scipy.fftpack import fft,ifft
import STA_LTA

class Quake():
    def __init__(self,path,data=None,delta=None,T_series=None):
        self.path = path
        self.data = read(self.path)[0].data
        self.delta = read(self.path)[0].stats.delta
        self.T_series = [i*self.delta for i in range(len(self.data))]

def fourier_trans(T,sr):
    complex_array = np.fft.fft(sr)
    y_ = np.fft.ifft(complex_array).real
    Freqency = np.fft.fftfreq(y_.size,T[1]-T[0])
    print(Freqency.T[262140:262150])
    Amplitude = np.abs(complex_array)
    Phase = np.angle(complex_array)
    return Freqency,Amplitude,Phase

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


filename = 'model1'
Path_0 = 'D:\\research\\dsm\\DSMsynTI-mpi-master\\models\\model1\\result1\\'+filename+'.Zs'
Path_1 = 'seis_data\\XB.ELYSE.00.HHR.D.2022.124.230000.SAC'
Path_set = [Path_0,Path_1]

def get_energy(n):
    Q = Quake(Path_set[n])
    data = Q.data
    T = Q.T_series
    delta = Q.delta
    Freqency,Amplitude,Phase = fourier_trans(T,data)
    nperiod = 60
    period_step = 0.2
    start_period = 10
    end_period = 40
    time_halfwidth = 1
    aut = 50.3
    nstep = round((end_period-start_period)/period_step)
    energy = []
    ylen = 1000
    
    for point in range(nstep):
        freqency, amplitude = filter(Freqency,Amplitude/2,start_period+point*period_step,time_halfwidth,aut)
        freq_signal = amplitude*np.cos(Phase) + 1j*amplitude*np.sin(Phase)
        filtered_signal = np.real(np.fft.ifft(freq_signal))
        envelope = obspy.signal.filter.envelope(filtered_signal)
        sampled_envelope = np.zeros((ylen))
        sample_step = floor(len(T)/ylen)
        for i in range(ylen):
            sampled_envelope[i] = envelope[sample_step*i]
        energy_of_env = [i**2 for i in sampled_envelope]
        energy.append(energy_of_env)

    plt.rcParams['xtick.direction'] = 'in' 
    plt.rcParams['ytick.direction'] = 'in' 
    plt.xlabel('Period(s)')
    plt.ylabel('Time after P(s)')
    plt.title('Energy(dB)')
    
    #translate to dB
    max_energy = np.max(energy)
    for i in range(nstep):
        for j in range(ylen):
            energy[i][j] = log10(energy[i][j]/max_energy)
    x = np.linspace(start_period,end_period,nstep)
    if n == 1:
        y = np.linspace(-1665.8,len(T)*delta-1665.8,ylen)
    elif n == 0:
        P_arrival_time = STA_LTA.get_P_arrival_time(Path_set[n])*delta
        y = np.linspace(-P_arrival_time,len(T)*delta-P_arrival_time,ylen)
    X,Y = np.meshgrid(x,y)
    Z = np.array(energy)
    Z = Z.T
    C=plt.contourf(X,Y,Z,levels = np.linspace(-10,0,100),cmap = 'jet')
    cbar=plt.colorbar()
    plt.clabel(C, inline=False, fontsize=10)
    plt.show()
    
def plot(path):
    signal = read(path)
    signal[0].plot()

def get_filtered_waveforms(n):
    Q = Quake(Path_set[n])
    data = Q.data
    T = Q.T_series
    delta = Q.delta
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
get_energy(0)
plot(Path_0)
