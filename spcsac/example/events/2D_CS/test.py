import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.colorbar as colorbar
import scipy.fft as ft
from matplotlib import colors
import numpy.fft as fft

vector = 'Z'
dic = '../event5/'
pyname = '2DCS_pos.py'
Original_file = dic+'Original_npy/' + '202R.npy'
pgv = np.load(Original_file)
M = np.size(pgv,0)
N = np.size(pgv,1)
im = ft.dct(ft.dct(pgv,type=2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
eps = 1e-7
