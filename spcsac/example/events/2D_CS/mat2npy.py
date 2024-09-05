import numpy as np
import scipy.io as sio

dic = '../event5/'
matfile = 'pos'

Posdic = sio.loadmat(dic+matfile+'.mat')
Position = np.array(Posdic[matfile])
I = np.ones(Position.shape)
np.save(dic+matfile+'.npy',Position-I)