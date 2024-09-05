from math import *
import numpy as np
import sys
import random
import scipy.fft as fft

np.seterr(divide='ignore',invalid='ignore')

Sampling_Rate = 0.08
Side_len = 5
fremax = 20

inputfile = sys.argv[1]
outputfile = sys.argv[2]
outputpoints = sys.argv[3]

Original_Signal = np.load(inputfile)
M = np.size(Original_Signal,0)
N = np.size(Original_Signal,1)
npoints = round(Sampling_Rate*Side_len**2)

def SquareJittered(M,N,Side_len,Original_Signal,npoints):
    rows = round(M/Side_len)
    cols = round(N/Side_len)
    allpoints = rows*cols*npoints
    Sampled_Signal = np.zeros((M,N))
    Position = np.zeros((allpoints,3))
    for i in range(rows):
        for j in range(cols):
            order = [k for k in range(Side_len**2)]
            random.shuffle(order)
            for k in range(npoints):
                select_point = order[k]
                row = i*Side_len + select_point//Side_len
                col = j*Side_len + select_point%Side_len
                num = npoints*i*cols + j*npoints + k
                Position[num][0] = row
                Position[num][1] = col
                Position[num][2] = Original_Signal[row][col]
    #Position = np.delete(Position,slice(0,100),0)
    return Position            
                
def CSOMP(ob,Original_Mat,t,Position,M,N):
    rk = ob
    err = np.linalg.norm(ob,ord=2)
    err = 0.01*err
    L = M*N
    x = np.zeros(L)
    A = np.zeros((t,L))
    An = np.zeros((t,L))
    Sk = np.zeros((t,t))
    Skpos = np.zeros(L,dtype=int)

    for j in range(L):
        x=j//N
        y=j%N
        mat = np.zeros((M,N))
        if x<=fremax and y<=fremax: 
            mat[x][y] = 1
        fre = fft.idct(fft.idct(mat,axis=-1,norm='ortho'),axis=-2,norm='ortho')
        for i in range(t):
            row = int(Position[i][0])
            col = int(Position[i][1])
            A[i][j] = fre[row][col]
    for col in range(L):
        An[:,col] = A[:,col]/np.linalg.norm(A[:,col])
    An[np.isnan(An)] = 0
    for k in range(L):
        res = np.dot(An.T,rk)
        maxindex = int(np.argmax(abs(res)))
        row = maxindex//N
        col = maxindex%N
        Skpos[k] = maxindex
        Sk[:,k] = A[:,Skpos[k]]
        Asm = np.dot(Sk[:,0:k+1].T,Sk[:,0:k+1])
        try:
            Asm = np.linalg.inv(Asm)
        except:
            #print(k,000)
            if k > t:
                continue
            else:
                return np.zeros((L))
        Asp = np.dot(Asm,Sk[:,0:k+1].T)
        rk = ob
        xk = np.dot(Asp,rk)
        rk = rk - np.dot(Sk[:,0:k+1],xk)
        ctn = np.linalg.norm(abs(rk),ord=2)
        An[:,Skpos[k]] = 0
        if ctn<err:
            #print(ctn/err)
            #print(k+1)
            xr = np.zeros(L)
            for i in range(k+1):
                xr[Skpos[i]] = xk[i]
            return xr
        if k > t:
            return np.zeros((L))
        
    
def MatMP(M,N,Position,Original):
    t = Position.shape[0]
    Onepiece = Position[:,2]
    Coe = CSOMP(Onepiece,Original,t,Position,M,N)
    Recovered_in_Sparse = np.zeros((M,N))
    for i in range(M):
        for j in range(N):
            try:
                Recovered_in_Sparse[i][j] = Coe[j+i*N]
            except:
                return np.zeros((M,N))
    return Recovered_in_Sparse

Recovered_in_Sparse = np.zeros((M,N))
while np.all(Recovered_in_Sparse==0):
    Position = SquareJittered(M,N,Side_len,Original_Signal,npoints)
    Recovered_in_Sparse = MatMP(M,N,Position,Original_Signal)
Recovered_Signal = fft.idct(fft.idct(Recovered_in_Sparse,axis=-1,norm='ortho'),axis=-2,norm='ortho')
np.save(outputfile,Recovered_Signal)
np.save(outputpoints,Position)

