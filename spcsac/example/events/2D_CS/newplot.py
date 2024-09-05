import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.colorbar as colorbar
import scipy.fft as ft
from matplotlib import colors
import numpy.fft as fft

dic = '../event5/'
Original_file1 = dic+'Original_npy/' + '200R.npy'
Original_file2 = dic+'Original_npy/' + '201R.npy'
Original_file3 = dic+'Original_npy/' + '202R.npy'
Recovered_file1 = dic+'Recovered_npy/'+ '200R_ts.npy'
Recovered_file2 = dic+'Recovered_npy/'+ '201R_ts.npy'
Recovered_file3 = dic+'Recovered_npy/'+ '202R_ts.npy'

pgv1 = np.load(Original_file1)
pgv2 = np.load(Original_file2)
pgv3 = np.load(Original_file3)
pgvrc1 = np.load(Recovered_file1)
pgvrc2 = np.load(Recovered_file2)
pgvrc3 = np.load(Recovered_file3)

f_t = 0
fig_nrow = 3
fig_ncol = 3
fig = plt.figure('Comparation',figsize=(10, 10))

ax1 = plt.subplot(fig_nrow,fig_ncol,1)
ax1.set_title("Original Wavefield", loc="center", fontsize=14) 
im1 = ft.dct(ft.dct(pgv1,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    image1 = plt.imshow(pgv1)
else:
    image1 = plt.imshow(im1)
ax2 = plt.subplot(fig_nrow,fig_ncol,2)
ax2.set_title("Recovered Wavefield", loc="center", fontsize=14) 
imrc1 = ft.dct(ft.dct(pgvrc1,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    imagerc1 = plt.imshow(pgvrc1)
else:
    imagerc1 = plt.imshow(imrc1)
ax3 = plt.subplot(fig_nrow,fig_ncol,3)
ax3.set_title("Residual Wavefield", loc="center", fontsize=14) 
mins1 = plt.imshow(abs(pgvrc1-pgv1)/pgv1,vmin=0,vmax=0.1)

ax4 = plt.subplot(fig_nrow,fig_ncol,4)
ax4.set_title("Original Wavefield", loc="center", fontsize=14) 
im1 = ft.dct(ft.dct(pgv1,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    image2 = plt.imshow(pgv2)
else:
    image2 = plt.imshow(im2)

ax5 = plt.subplot(fig_nrow,fig_ncol,5)
ax5.set_title("Recovered Wavefield", loc="center", fontsize=14) 
imrc2 = ft.dct(ft.dct(pgvrc2,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    imagerc2 = plt.imshow(pgvrc2)
else:
    imagerc2 = plt.imshow(imrc2)

ax6 = plt.subplot(fig_nrow,fig_ncol,6)
ax6.set_title("Residual Wavefield", loc="center", fontsize=14) 
mins2 = plt.imshow(abs(pgvrc2-pgv2)/pgv2)

ax7 = plt.subplot(fig_nrow,fig_ncol,7)
ax7.set_title("Original Wavefield", loc="center", fontsize=14) 
im1 = ft.dct(ft.dct(pgv1,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    image3 = plt.imshow(pgv3)
else:
    image3 = plt.imshow(im3)

ax8 = plt.subplot(fig_nrow,fig_ncol,8)
ax8.set_title("Recovered Wavefield", loc="center", fontsize=14) 
imrc3 = ft.dct(ft.dct(pgvrc3,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
if f_t==0:
    imagerc3 = plt.imshow(pgvrc3)
else:
    imagerc3 = plt.imshow(imrc3)
    
ax9 = plt.subplot(fig_nrow,fig_ncol,9)
ax9.set_title("Residual Wavefield", loc="center", fontsize=14) 
mins3 = plt.imshow(abs(pgvrc3-pgv3)/pgv3)
plt.show()