import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.colorbar as colorbar
import scipy.fft as ft
from matplotlib import colors
import subprocess
import numpy.fft as fft
import noise

vector = 'Z'
dic = '../event5/'
pyname = '2DCS_pos.py'
Original_file = dic+'Original_npy/' + '202R.npy'
pgv = np.load(Original_file)
# #pgv = noise.Noise(pgv)
# Noised_file = dic + 'Original_npy/'+ vector +'Noise_pgv.npy'
# np.save(Noised_file,pgv)
Recovered_file = dic+'Recovered_npy/'+ vector +'pgvrc.npy'
Position_file = dic+'pos.npy'
position = np.load(Position_file)
x = position[:,1]
y = position[:,0]
x[-1] = 57
y[-1] = 6
x[-2] = 55
y[-2] = 40
np.save(Position_file,position)
info_file = dic + 'info.txt'
#ev5:184,300
srow = '184'
scol = '300'
subprocess.run(['python',pyname,Original_file,Recovered_file,Position_file,srow,scol])


info = np.loadtxt(info_file,skiprows=0,dtype='float',comments='#')
source_lon,source_lat,lon_s,lon_e,lat_s,lat_e,X_num,Y_num = info
xtick_num = 7
ytick_num = 6
lon_step = (lon_e-lon_s)/(xtick_num-1)
lat_step = (lat_e-lat_s)/(ytick_num-1)
X_step = X_num/(xtick_num-1)
Y_step = Y_num/(ytick_num-1)
pgvrc = np.load(Recovered_file)
vmin = min(np.min(pgv), np.min(pgvrc))
vmax = max(np.max(pgv), np.max(pgvrc))
norm = colors.Normalize(vmin=vmin, vmax=vmax)

f_t = 0
show_station = 1
fig_nrow = 1
fig_ncol = 3
cmap = 'Spectral_r'
fig = plt.figure('Comparation',figsize=(10, 10))
fig.suptitle("{} Component   Sampling Rate = 0.064".format(vector),x=0.5,y=0.8, fontsize=16)
ax1 = plt.subplot(fig_nrow,fig_ncol,1)
ax1.set_title("Original Wavefield", loc="center", fontsize=14) 
im = ft.dct(ft.dct(pgv,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    image = plt.imshow(pgv,vmin = vmin, vmax = vmax,cmap=cmap)
else:
    image = plt.imshow(im,vmin = vmin, vmax = vmax,cmap=cmap)
if show_station == 1:
    s1 = plt.scatter(x,y,s=10, c='k', marker='^')
#s2 = plt.scatter(16.67,16.67,s=100, c='red', marker='*')
#plt.legend((s1,s2),('Sampled Station','epicenter'),bbox_to_anchor=(0.2, 0.8),bbox_transform=plt.gcf().transFigure)
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
ax2 = plt.subplot(fig_nrow,fig_ncol,2)
ax2.set_title("Recovered Wavefield", loc="center", fontsize=14) 
im2 = ft.dct(ft.dct(pgvrc,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    imagerc = plt.imshow(pgvrc,vmin = vmin, vmax = vmax,cmap=cmap)
else:
    imagerc = plt.imshow(im2,vmin = vmin, vmax = vmax,cmap=cmap)
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
plt.scatter(16.67,16.67,s=100, c='red', marker='*')
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
cax = fig.add_axes([0.126, 0.205, 0.501, 0.02])
plt.colorbar(imagerc,cax=cax,label='Velocity(m/s)',orientation="horizontal")

ax3 = plt.subplot(fig_nrow,fig_ncol,3)
ax3.set_title("Residual Wavefield", loc="center", fontsize=14) 
mins = plt.imshow(abs(pgvrc-pgv)/pgv,vmin=0,vmax=0.1,cmap=cmap)
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
plt.colorbar(mins,label='err',orientation="horizontal")

err = 0
value = 0
for i in range(len(x)):
    row = int(position[i][0])
    col = int(position[i][1])
    err = err + ((pgv[row][col] - pgvrc[row][col]))**2
    value = value + pgv[row][col]**2
value = (value/len(x))**0.5
err = (err/len(x))**0.5/value
err2 = np.sum(np.square(pgvrc-pgv))
value2 = np.sum(np.square(pgv))
err2 = (err2/value2)**0.5
print(err,err2)
plt.show()
