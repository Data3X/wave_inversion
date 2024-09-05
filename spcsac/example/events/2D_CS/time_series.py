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
pyname = '2DCS_BB.py'
Original_file1 = dic+'Original_npy/' + '200R.npy'
Original_file2 = dic+'Original_npy/' + '201R.npy'
Original_file3 = dic+'Original_npy/' + '202R.npy'
pgv1 = np.load(Original_file1)
pgv2 = np.load(Original_file2)
pgv3 = np.load(Original_file3)
# #pgv = noise.Noise(pgv)
# Noised_file = dic + 'Original_npy/'+ vector +'Noise_pgv.npy'
# np.save(Noised_file,pgv)
Recovered_file1 = dic+'Recovered_npy/'+ '200R_ts.npy'
Recovered_file2 = dic+'Recovered_npy/'+ '201R_ts.npy'
Recovered_file3 = dic+'Recovered_npy/'+ '202R_ts.npy'
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
subprocess.run(['python',pyname,Original_file1,Original_file2,Original_file3,Recovered_file1,Recovered_file2,Recovered_file3,Position_file])


info = np.loadtxt(info_file,skiprows=0,dtype='float',comments='#')
source_lon,source_lat,lon_s,lon_e,lat_s,lat_e,X_num,Y_num = info
xtick_num = 7
ytick_num = 6
lon_step = (lon_e-lon_s)/(xtick_num-1)
lat_step = (lat_e-lat_s)/(ytick_num-1)
X_step = X_num/(xtick_num-1)
Y_step = Y_num/(ytick_num-1)
pgvrc1 = np.load(Recovered_file1)
pgvrc2 = np.load(Recovered_file2)
pgvrc3 = np.load(Recovered_file3)
vmin = min(np.min(pgv1), np.min(pgvrc1),np.min(pgv2), np.min(pgvrc2),np.min(pgv3), np.min(pgvrc3))
vmax = max(np.max(pgv1), np.max(pgvrc1),np.max(pgv2), np.max(pgvrc2),np.max(pgv3), np.max(pgvrc3))
norm = colors.Normalize(vmin=vmin, vmax=vmax)

f_t = 0
show_station = 1
fig_nrow = 3
fig_ncol = 3
cmap = 'Spectral_r'
fig = plt.figure('Comparation',figsize=(10, 10))
fig.suptitle("{} Component   Sampling Rate = 0.064".format(vector),x=0.5,y=0.8, fontsize=16)
ax1 = plt.subplot(fig_nrow,fig_ncol,1)
ax1.set_title("Original Wavefield", loc="center", fontsize=14) 
im1 = ft.dct(ft.dct(pgv1,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    image1 = plt.imshow(pgv1,vmin = vmin, vmax = vmax,cmap=cmap)
else:
    image1 = plt.imshow(im1,vmin = vmin, vmax = vmax,cmap=cmap)
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
#s2 = plt.scatter(16.67,16.67,s=100, c='red', marker='*')
#plt.legend((s1,s2),('Sampled Station','epicenter'),bbox_to_anchor=(0.2, 0.8),bbox_transform=plt.gcf().transFigure)
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
ax2 = plt.subplot(fig_nrow,fig_ncol,2)
ax2.set_title("Recovered Wavefield", loc="center", fontsize=14) 
imrc1 = ft.dct(ft.dct(pgvrc1,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    imagerc1 = plt.imshow(pgvrc1,vmin = vmin, vmax = vmax,cmap=cmap)
else:
    imagerc1 = plt.imshow(imrc1,vmin = vmin, vmax = vmax,cmap=cmap)
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
plt.scatter(16.67,16.67,s=100, c='red', marker='*')
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
cax = fig.add_axes([0.126, 0.205, 0.501, 0.02])
plt.colorbar(imagerc1,cax=cax,label='Velocity(m/s)',orientation="horizontal")

ax3 = plt.subplot(fig_nrow,fig_ncol,3)
ax3.set_title("Residual Wavefield", loc="center", fontsize=14) 
mins1 = plt.imshow(abs(pgvrc1-pgv1)/pgv1,vmin=0,vmax=0.1,cmap=cmap)
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
plt.colorbar(mins1,label='err',orientation="horizontal")
############################################################################
ax4 = plt.subplot(fig_nrow,fig_ncol,4)
ax4.set_title("Original Wavefield", loc="center", fontsize=14) 
im1 = ft.dct(ft.dct(pgv1,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    image2 = plt.imshow(pgv2,vmin = vmin, vmax = vmax,cmap=cmap)
else:
    image2 = plt.imshow(im2,vmin = vmin, vmax = vmax,cmap=cmap)
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
#s2 = plt.scatter(16.67,16.67,s=100, c='red', marker='*')
#plt.legend((s1,s2),('Sampled Station','epicenter'),bbox_to_anchor=(0.2, 0.8),bbox_transform=plt.gcf().transFigure)
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])

ax5 = plt.subplot(fig_nrow,fig_ncol,5)
ax5.set_title("Recovered Wavefield", loc="center", fontsize=14) 
imrc2 = ft.dct(ft.dct(pgvrc2,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    imagerc2 = plt.imshow(pgvrc2,vmin = vmin, vmax = vmax,cmap=cmap)
else:
    imagerc2 = plt.imshow(imrc2,vmin = vmin, vmax = vmax,cmap=cmap)
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
plt.scatter(16.67,16.67,s=100, c='red', marker='*')
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
cax = fig.add_axes([0.126, 0.205, 0.501, 0.02])
plt.colorbar(imagerc2,cax=cax,label='Velocity(m/s)',orientation="horizontal")

ax6 = plt.subplot(fig_nrow,fig_ncol,6)
ax6.set_title("Residual Wavefield", loc="center", fontsize=14) 
mins2 = plt.imshow(abs(pgvrc2-pgv2)/pgv2,vmin=0,vmax=0.1,cmap=cmap)
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
plt.colorbar(mins2,label='err',orientation="horizontal")
###################################################################################
ax7 = plt.subplot(fig_nrow,fig_ncol,7)
ax7.set_title("Original Wavefield", loc="center", fontsize=14) 
im1 = ft.dct(ft.dct(pgv1,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    image3 = plt.imshow(pgv3,vmin = vmin, vmax = vmax,cmap=cmap)
else:
    image3 = plt.imshow(im3,vmin = vmin, vmax = vmax,cmap=cmap)
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
#s2 = plt.scatter(16.67,16.67,s=100, c='red', marker='*')
#plt.legend((s1,s2),('Sampled Station','epicenter'),bbox_to_anchor=(0.2, 0.8),bbox_transform=plt.gcf().transFigure)
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
ax8 = plt.subplot(fig_nrow,fig_ncol,8)
ax8.set_title("Recovered Wavefield", loc="center", fontsize=14) 
imrc3 = ft.dct(ft.dct(pgvrc3,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    imagerc3 = plt.imshow(pgvrc3,vmin = vmin, vmax = vmax,cmap=cmap)
else:
    imagerc3 = plt.imshow(imrc3,vmin = vmin, vmax = vmax,cmap=cmap)
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
plt.scatter(16.67,16.67,s=100, c='red', marker='*')
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
cax = fig.add_axes([0.126, 0.205, 0.501, 0.02])
plt.colorbar(imagerc3,cax=cax,label='Velocity(m/s)',orientation="horizontal")

ax9 = plt.subplot(fig_nrow,fig_ncol,9)
ax9.set_title("Residual Wavefield", loc="center", fontsize=14) 
mins3 = plt.imshow(abs(pgvrc3-pgv3)/pgv3,vmin=0,vmax=0.1,cmap=cmap)
plt.xticks(np.arange(0,X_num,X_step),['{}\u00B0E'.format(round(i,2)) for i in np.arange(lon_e,lon_s,-lon_step)])
plt.yticks(np.arange(0,Y_num,Y_step),['{}\u00B0N'.format(round(i,2)) for i in np.arange(lat_e,lat_s,-lat_step)])
if show_station == 1:
    plt.scatter(x,y,s=10, c='k', marker='^')
plt.colorbar(mins3,label='err',orientation="horizontal")
plt.show()
