from obspy import read
from matplotlib import pyplot as plt
from math import *
import scipy

def get_P_arrival_time(filepath):
    singlechannel = read(filepath)
    data = singlechannel[0].data
    length = len(data)
    long_window = 200
    short_window = 20
    sta = 0
    lta = 0
    STA_LTA = []
    def func(start,t):
        return (data[t+start+1])**2

    for start in range(length-long_window-short_window): 
        for lcount in range(long_window):
            lta = lta + func(start,lcount)
        for scount in range(short_window):
            sta = sta + func(start+long_window,scount)
        STA_LTA.append(sta*long_window/lta/short_window-1)

    plt.figure(1)
    plt.plot(STA_LTA)
    plt.show()

    amp = max(STA_LTA)
    if amp > 1:
        time = long_window+short_window/2+STA_LTA.index(amp)
        return time

    singlechannel.plot()
filename = 'model1'
Path_0 = 'D:\\research\\dsm\\DSMsynTI-mpi-master\\models\\model1\\result1\\'+filename+'.Zs'
read(Path_0)[0].plot()
