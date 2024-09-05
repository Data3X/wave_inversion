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

def CSOMP(ob1,ob2,ob3,Position,M,N,t,srow=None,scol=None):
    lambda1 = 0.0
    lambda2 = 0.0
    ob = np.concatenate((ob1,ob2,ob3),axis=0)
    rk1 = ob1
    rk2 = ob2
    rk3 = ob3
    err1 = np.linalg.norm(ob1,ord=2)
    err2 = np.linalg.norm(ob2,ord=2)
    err3 = np.linalg.norm(ob3,ord=2)
    err1 = 0.01*err1
    err2 = 0.01*err2
    err3 = 0.01*err3
    L = M*N
    kmax = 3*t
    x = np.zeros(L)
    A = np.zeros((5*t,L))
    B1 = np.zeros((t,L))
    B2 = np.zeros((t,L))
    B3 = np.zeros((t,L))
    An = np.zeros((t,L))
    Sk1 = np.zeros((t,kmax))
    Sk2 = np.zeros((t,kmax))
    Sk3 = np.zeros((t,kmax))
    Sk2S2 = np.zeros((t,kmax))
    Sk2S3 = np.zeros((t,kmax))
    Sk2S4 = np.zeros((t,kmax))
    Sk2S5 = np.zeros((t,kmax))
    Skpos1 = np.zeros(L,dtype=int)
    Skpos2 = np.zeros(L,dtype=int)
    Skpos3 = np.zeros(L,dtype=int)

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
            A[i][j] = u[row][col]
            A[i+t][j] = u[row][col+1]
            A[i+2*t][j] = u[row-1][col]
            A[i+3*t][j] = u[row][col-1]
            A[i+4*t][j] = u[row+1][col]
        An[:,j] = A[0:t,j]/np.linalg.norm(A[0:t,j])
        fre_weight = min((1-x/M),(1-y/N))
        for i in range(t): 
            row = int(Position[i][0])
            col = int(Position[i][1])
            B1[i][j] = An[i][j]*dirmat[row][col]*fre_weight
    B2 = B1
    B3 = B1
    B1[np.isnan(An)] = 0
    B2[np.isnan(An)] = 0
    B3[np.isnan(An)] = 0
    for k in range(L):
        v1 = np.dot(B1.T,rk1)
        v2 = np.dot(B2.T,rk2)
        v3 = np.dot(B3.T,rk3)
        maxindex1 = int(np.argmax(abs(v1)))
        maxindex2 = int(np.argmax(abs(v2)))
        maxindex3 = int(np.argmax(abs(v3)))
        row1 = maxindex1//N
        row2 = maxindex2//N
        row3 = maxindex3//N
        col1 = maxindex1%N
        col2 = maxindex2%N
        col3 = maxindex3%N
        Skpos1[k] = maxindex1
        Skpos2[k] = maxindex2
        Skpos3[k] = maxindex3
        Sk1[:,k] = A[0:t,Skpos1[k]]
        Sk2[:,k] = A[0:t,Skpos2[k]]
        Sk3[:,k] = A[0:t,Skpos3[k]]
        Sk2S2[:,k] = A[t:2*t,Skpos2[k]]
        Sk2S3[:,k] = A[2*t:3*t,Skpos2[k]]
        Sk2S4[:,k] = A[3*t:4*t,Skpos2[k]]
        Sk2S5[:,k] = A[4*t:5*t,Skpos2[k]]
        d1 = np.dot(Sk1[:,0:k+1].T,ob1)
        d2 = np.dot(Sk1[:,0:k+1].T,ob2)
        d3 = np.dot(Sk1[:,0:k+1].T,ob3)
        d = np.concatenate((d1,d2,d3),axis=0)
        C = np.zeros((3*(k+1),3*(k+1)))
        C[0:k+1,0:k+1] = np.dot(Sk1[:,0:k+1].T,Sk1[:,0:k+1])
        C[k+1:2*k+2,k+1:2*k+2] = np.dot(Sk2[:,0:k+1].T,Sk2[:,0:k+1])
        C[2*k+2:3*k+3,2*k+2:3*k+3] = np.dot(Sk3[:,0:k+1].T,Sk3[:,0:k+1])
        mat1 = lambda1*np.concatenate((Sk1[:,0:k+1].T,-Sk2S2[:,0:k+1].T,np.zeros((k+1,t))),axis=0)
        mat2 = np.concatenate((Sk1[:,0:k+1],-Sk2S2[:,0:k+1],np.zeros((t,k+1))),axis=1)
        mat3 = lambda2*np.concatenate((np.zeros((k+1,t)),-Sk2S2[:,0:k+1].T,Sk3[:,0:k+1].T),axis=0)
        mat4 = np.concatenate((np.zeros((t,k+1)),-Sk2S2[:,0:k+1],Sk3[:,0:k+1]),axis=1)
        C = C + np.dot(mat1,mat2) + np.dot(mat3,mat4)
        iC = np.linalg.inv(C)
        xk = np.dot(iC,d)
        rk1 = np.dot(Sk1[:,0:k+1],xk[0:k+1,0])
        rk2 = np.dot(Sk2[:,0:k+1],xk[k+1:2*k+2,0])
        rk3 = np.dot(Sk3[:,0:k+1],xk[2*k+2:3*k+3,0])
        rg12 = np.dot(Sk2S2[:,0:k+1],xk[0:k+1,0])
        rg23 = np.dot(Sk2S2[:,0:k+1],xk[0:k+1,0])
        res1 = rk1 - ob1[:,0]
        res2 = rk2 - ob2[:,0]
        res3 = rk3 - ob3[:,0]
        res4 = lambda1*(rk1 - rg12)
        res5 = lambda2*(rk3 - rg23)
        res = np.concatenate((res1,res2,res3,res4,res5),axis=0)
        ctn = np.linalg.norm(abs(res),ord=2)/np.linalg.norm(abs(ob),ord=2)
        B1[:,Skpos1[k]] = 0
        B2[:,Skpos2[k]] = 0
        B3[:,Skpos3[k]] = 0
        if ctn<0.001 or k>=2*t:
            xr1 = np.zeros(L)
            xr2 = np.zeros(L)
            xr3 = np.zeros(L)
            for i in range(k):
                xr1[Skpos1[i]] = xk[i]
                xr2[Skpos2[i]] = xk[i+k+1]
                xr3[Skpos3[i]] = xk[i+2*k+2]
            return xr1,xr2,xr3
        
def MatMP(M,N,Position,Original1,Original2,Original3,srow=None,scol=None):
    t = Position.shape[0]
    Onepiece1 = np.zeros((t,1))
    Onepiece2 = np.zeros((t,1))
    Onepiece3 = np.zeros((t,1))
    for i in range(t):
        row = int(Position[i][0])
        col = int(Position[i][1])
        Onepiece1[i][0] = Original1[row][col]
        Onepiece2[i][0] = Original2[row][col]
        Onepiece3[i][0] = Original3[row][col]
    Coe1,Coe2,Coe3 = CSOMP(Onepiece1,Onepiece2,Onepiece3,Position,M,N,t,srow,scol)
    Recovered_in_Sparse1,Recovered_in_Sparse2,Recovered_in_Sparse3 = np.zeros((M,N)),np.zeros((M,N)),np.zeros((M,N))
    for i in range(M):
        for j in range(N):
            try:
                Recovered_in_Sparse1[i][j] = Coe1[j+i*N]
                Recovered_in_Sparse2[i][j] = Coe2[j+i*N]
                Recovered_in_Sparse3[i][j] = Coe3[j+i*N]
            except:
                return np.zeros((M,N)),np.zeros((M,N)),np.zeros((M,N))
    return Recovered_in_Sparse1,Recovered_in_Sparse2,Recovered_in_Sparse3

np.seterr(divide='ignore',invalid='ignore')
xfremax = 23
yfremax = 30
inputfile1 = sys.argv[1]
inputfile2 = sys.argv[2]
inputfile3 = sys.argv[3]
outputfile1 = sys.argv[4]
outputfile2 = sys.argv[5]
outputfile3 = sys.argv[6]
pos = sys.argv[7]
try:
    srow = int(sys.argv[8])
    scol = int(sys.argv[9])
    r = int(sys.argv[10])
except:
    srow = None
    scol = None
    r = None
Original_Signal1 = np.load(inputfile1,allow_pickle=True)
Original_Signal2 = np.load(inputfile2,allow_pickle=True)
Original_Signal3 = np.load(inputfile3,allow_pickle=True)
M = np.size(Original_Signal1,0)
N = np.size(Original_Signal1,1)   
Position = np.load(pos,allow_pickle=True)
Recovered_in_Sparse1,Recovered_in_Sparse2,Recovered_in_Sparse3 = MatMP(M,N,Position,Original_Signal1,Original_Signal2,Original_Signal3,srow,scol)
Recovered_Signal1 = fft.idct(fft.idct(Recovered_in_Sparse1,axis=-1,norm='ortho'),axis=-2,norm='ortho')
Recovered_Signal2 = fft.idct(fft.idct(Recovered_in_Sparse2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
Recovered_Signal3 = fft.idct(fft.idct(Recovered_in_Sparse3,axis=-1,norm='ortho'),axis=-2,norm='ortho')
np.save(outputfile1,Recovered_Signal1)
np.save(outputfile2,Recovered_Signal2)
np.save(outputfile3,Recovered_Signal3)

