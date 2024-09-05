import obspy
from math import *
import numpy as np
import matplotlib.pyplot as plt

def getlog(x):
    if x>0:
        return log10(x)
    elif x<0:
        return log10(-x)
    else:
        return x

nrow = 46
ncol = 61
num = nrow*ncol
date = '.200710310330A'
dic = '../event3/'
Rtrace = [[[] for i in range(ncol)] for i in range(nrow)]
Ttrace = [[[] for i in range(ncol)] for i in range(nrow)]
Ztrace = [[[] for i in range(ncol)] for i in range(nrow)]
for i in range(nrow):
    for j in range(ncol):
        num = i*ncol+j+1
        Rf = dic + 'sacfile/' + str(num)+date+'.Rs'
        Tf = dic + 'sacfile/' + str(num)+date+'.Ts'
        Zf = dic + 'sacfile/' + str(num)+date+'.Zs'
        Rtrace[i][j] = obspy.read(Rf)[0].data.tolist()
        Ttrace[i][j] = obspy.read(Tf)[0].data.tolist()
        Ztrace[i][j] = obspy.read(Zf)[0].data.tolist()

Rarray = np.array(Rtrace)
Tarray = np.array(Ttrace)
Zarray = np.array(Ztrace)
length = len(Rtrace[0][0])
RTarray = np.zeros((nrow,ncol,length))
RZarray = np.zeros((nrow,ncol,length))
RTZarray = np.zeros((nrow,ncol,length))
Rpgv = np.zeros((nrow,ncol))
Tpgv = np.zeros((nrow,ncol))
Zpgv = np.zeros((nrow,ncol))
RTpgv = np.zeros((nrow,ncol))
RZpgv = np.zeros((nrow,ncol))
RTZpgv = np.zeros((nrow,ncol))


for i in range(nrow):
    for j in range(ncol):
        Rmax = 0
        Tmax = 0
        Zmax = 0
        RTmax = 0
        RZmax = 0
        RTZmax = 0
        for k in range(length):
            RTarray[i][j][k] = sqrt(Rarray[i][j][k]**2+Tarray[i][j][k]**2)
            RZarray[i][j][k] = sqrt(Rarray[i][j][k]**2+Zarray[i][j][k]**2)
            RTZarray[i][j][k] = sqrt(Rarray[i][j][k]**2+Tarray[i][j][k]**2+Zarray[i][j][k]**2)
            if Rarray[i][j][k] > Rmax:
                Rmax = Rarray[i][j][k]
            if Tarray[i][j][k] > Tmax:
                Tmax = Tarray[i][j][k]
            if Zarray[i][j][k] > Zmax:
                Zmax = Zarray[i][j][k]
            if RTarray[i][j][k] > RTmax:
                RTmax = RTarray[i][j][k]
            if RZarray[i][j][k] > RZmax:
                RZmax = RZarray[i][j][k]
            if RTZarray[i][j][k] > RTZmax:
                RTZmax = RTZarray[i][j][k]
            # Rarray[i][j][k] = getlog(Rarray[i][j][k])
            # Tarray[i][j][k] = getlog(Tarray[i][j][k])
            # Zarray[i][j][k] = getlog(Zarray[i][j][k])
        Rpgv[i][j] = Rmax
        Tpgv[i][j] = Tmax
        Zpgv[i][j] = Zmax
        RTpgv[i][j] = RTmax
        RZpgv[i][j] = RZmax
        RTZpgv[i][j] = RTZmax
    
suffix = 'npy'
for i in range(length):
    num = str(i+1)
    R = Rarray[:,:,i]
    T = Tarray[:,:,i]
    Z = Zarray[:,:,i]
    RT = RTarray[:,:,i]
    RZ = RZarray[:,:,i]
    RTZ = RTZarray[:,:,i]
    np.save(dic + 'Original_npy/' + num+'R.'+suffix,R)
    np.save(dic + 'Original_npy/' + num+'T.'+suffix,T)
    np.save(dic + 'Original_npy/' + num+'Z.'+suffix,Z)
    np.save(dic + 'Original_npy/' + num+'RT.'+suffix,RT)
    np.save(dic + 'Original_npy/' + num+'RZ.'+suffix,RZ)
    np.save(dic + 'Original_npy/' + num+'RTZ.'+suffix,RTZ)

np.save(dic + 'Original_npy/' + 'Rpgv.'+suffix,Rpgv)
np.save(dic + 'Original_npy/' + 'Tpgv.'+suffix,Tpgv)
np.save(dic + 'Original_npy/' + 'Zpgv.'+suffix,Zpgv)
np.save(dic + 'Original_npy/' + 'RTpgv.'+suffix,RTpgv)
np.save(dic + 'Original_npy/' + 'RZpgv.'+suffix,RZpgv)
np.save(dic + 'Original_npy/' + 'RTZpgv.'+suffix,RTZpgv)

plt.imshow(Zpgv)
plt.show()