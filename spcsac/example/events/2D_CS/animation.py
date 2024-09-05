import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.colorbar as colorbar
import scipy.fft as ft
from matplotlib import colors
import scipy.io as sio

Component = 'T'
dic = '../event3/'
Position_file = dic+'pos.npy'
info_file = dic + 'info.txt'

position = np.load(Position_file)
x = position[:,1]
y = position[:,0]

info = np.loadtxt(info_file,skiprows=0,dtype='float',comments='#')
source_lon,source_lat,lon_s,lon_e,lat_s,lat_e,X_num,Y_num = info
xtick_num = 7
ytick_num = 6
lon_step = (lon_e-lon_s)/(xtick_num-1)
lat_step = (lat_e-lat_s)/(ytick_num-1)
X_step = X_num/(xtick_num-1)
Y_step = Y_num/(ytick_num-1)

frame_num = 2048
delta = 1/16
period = round(frame_num*delta,2)

fig = plt.figure('Comparation',figsize=(10, 10))
for i in range(1,frame_num):
    time = round((i+1)*delta,2)
    fig.suptitle("Wavefield of {} Component   Time:{}/{}s".format(Component,time,period),x=0.5,y=0.8, fontsize=16)
    frameid = dic + 'Original_npy/' + str(i+1)+ Component + '.npy'
    or_field = np.load(frameid)
    reframeid =dic + 'Recovered_npy/' + str(i+1)+'Recovered_' + Component + '.npy'
    try:
        re_field = np.load(reframeid)
    except:
        pass
    # vmax = 1e-7
    # vmin = -vmax
    # norm = colors.Normalize(vmin=vmin, vmax=vmax)
    ######################
    ax1 = plt.subplot(131)
    ax1.set_title("Original Wavefield", loc="center", fontsize=14) 
    image = plt.imshow(or_field)
    plt.scatter(x,y,s=10, c='k', marker='^')
    plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_s,lon_e,lon_step)])
    plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_s,lat_e,lat_step)])
    ######################
    ax2 = plt.subplot(132)
    ax2.set_title("Recovered Wavefield", loc="center", fontsize=14) 
    try:
        reimage = plt.imshow(re_field)
    except:
        pass
    s1 = plt.scatter(x,y,s=10, c='k', marker='^')
    plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_s,lon_e,lon_step)])
    plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_s,lat_e,lat_step)])
    ######################
    ax3 = plt.subplot(133)
    ax3.set_title("Residual Wavefield", loc="center", fontsize=14) 
    try:
        mins = plt.imshow(abs(re_field-or_field),vmin = vmin, vmax = vmax)
    except:
        pass
    plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_s,lon_e,lon_step)])
    plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_s,lat_e,lat_step)])
    a = ax1.pcolormesh(or_field, cmap=plt.get_cmap('rainbow'))
    try:
        b = ax2.pcolormesh(re_field,norm=norm,cmap=plt.get_cmap('rainbow'))
    except:
        pass
    cax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
    fig.colorbar(a,cax=cax,ax=[ax1,ax2,ax3],label='Velocity(m/s)')
    plt.pause(0.01)
    plt.clf()
