from math import *
import numpy as np
import sys
import random
import scipy.fft as fft
import scipy.io as sio     

def getdirection(Position,t,M,N,i,j,srow,scol,r=None):
    if i==0:
        alphap = 1/sqrt(M)
    else:
        alphap = sqrt(2/M)
    if j==0:
        alphaq = 1/sqrt(N)
    else:
        alphaq = sqrt(2/N)
    mat = np.ones((M,N))
    weight_type = 1
    for k in range(t):
        m = int(Position[k][0])
        n = int(Position[k][1])
        pApm = -alphap*alphaq*(pi*i/M)*sin(pi*(2*m+1)*i/(2*M))*cos(pi*(2*n+1)*j/(2*N))
        pApn = -alphap*alphaq*(pi*j/N)*sin(pi*(2*n+1)*j/(2*N))*cos(pi*(2*m+1)*i/(2*M))
        if (scol-n)**2+(srow-m)**2 > r**2:
            if pApm == 0 and pApn==0:
                mat[m][n] = 1
            else: 
                a = np.array([pApm,pApn])
                b = np.array([srow-m,scol-n])
                if weight_type==0:
                    mat[m][n] = abs(np.dot(a,b))/(np.linalg.norm(a)*np.linalg.norm(b))
                else:
                    mat[m][n] = 1-2*acos(abs(np.dot(a,b))/(np.linalg.norm(a)*np.linalg.norm(b)))/pi
    return mat

def CSOMP(ob,Position,M,N,t,srow=None,scol=None):
    rk = ob
    err = np.linalg.norm(ob,ord=2)
    err = 0.01*err
    L = M*N
    x = np.zeros(L)
    A = np.zeros((t,L))
    B = np.zeros((t,L))
    An = np.zeros((t,L))
    Sk = np.zeros((t,t))
    Skpos = np.zeros(L,dtype=int)

    for j in range(L):
        x=j//N
        y=j%N
        mat = np.zeros((M,N))
        dirmat = np.ones((M,N))
        if x<=xfremax and y<=yfremax: 
            mat[x][y] = 1
            if srow:
                dirmat = getdirection(Position,t,M,N,x,y,srow,scol,4)
        u = fft.idct(fft.idct(mat,type=2,axis=-1,norm='ortho'),type=2,axis=-2,norm='ortho')
        for i in range(t):
            row = int(Position[i][0])
            col = int(Position[i][1])
            #print(dirmat[row][col])
            A[i][j] = u[row][col]
        An[:,j] = A[:,j]/np.linalg.norm(A[:,j])
        fre_weight = min((1-x/M),(1-y/N))
        for i in range(t):
            row = int(Position[i][0])
            col = int(Position[i][1])
            B[i][j] = An[i][j]*dirmat[row][col]*fre_weight

    for k in range(L):
        res = np.dot(B.T,rk)
        maxindex = int(np.argmax(abs(res)))
        row = maxindex//N
        col = maxindex%N
        Skpos[k] = maxindex
        Sk[:,k] = A[:,Skpos[k]]
        try:
            Asm = np.dot(Sk[:,0:k+1].T,Sk[:,0:k+1])
            Asm = np.linalg.inv(Asm)
        except:
            Asm = np.dot(Sk[:,0:k+1].T,Sk[:,0:k+1]) 
            lam = np.sum(abs(Asm))/M/N/10
            Asm = Asm + lam*np.eye(k+1)
            Asm = np.linalg.inv(Asm)
        Asp = np.dot(Asm,Sk[:,0:k+1].T)
        rk = ob
        xk = np.dot(Asp,rk)
        rk = rk - np.dot(Sk[:,0:k+1],xk)
        ctn = np.linalg.norm(abs(rk),ord=2)
        B[:,Skpos[k]] = 0
        if ctn<err or k>=8:
            #print(ctn/err)
            #print(k+1)
            xr = np.zeros(L)
            for i in range(k+1):
                xr[Skpos[i]] = xk[i]
            return xr
        
def MatMP(M,N,Position,Original,srow=None,scol=None):
    t = Position.shape[0]
    Onepiece = np.zeros(t)
    for i in range(t):
        row = int(Position[i][0])
        col = int(Position[i][1])
        Onepiece[i] = Original[row][col]
    Coe = CSOMP(Onepiece,Position,M,N,t,srow,scol)
    Recovered_in_Sparse = np.zeros((M,N))
    for i in range(M):
        for j in range(N):
            try:
                Recovered_in_Sparse[i][j] = Coe[j+i*N]
            except:
                return np.zeros((M,N))
    return Recovered_in_Sparse

np.seterr(divide='ignore',invalid='ignore')
xfremax = 23
yfremax = 30
inputfile = sys.argv[1]
outputfile = sys.argv[2]
pos = sys.argv[3]
try:
    srow = int(sys.argv[4])
    scol = int(sys.argv[5])
    r = int(sys.argv[6])
except:
    srow = None
    scol = None
    r = None
Original_Signal = np.load(inputfile,allow_pickle=True)
M = np.size(Original_Signal,0)
N = np.size(Original_Signal,1)   
Position = np.load(pos,allow_pickle=True)
Recovered_in_Sparse = MatMP(M,N,Position,Original_Signal,srow,scol)
Recovered_Signal = fft.idct(fft.idct(Recovered_in_Sparse,axis=-1,norm='ortho'),axis=-2,norm='ortho')
np.save(outputfile,Recovered_Signal)

