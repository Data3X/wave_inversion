from math import *
import numpy as np
import sys
import random
import scipy.fft as fft
import scipy.io as sio     
import matplotlib.pyplot as plt

def getOriginalSamplingPointsValue(OS):
    ob = np.zeros(t)
    for i in range(t):
        row = int(Position[i][0])
        col = int(Position[i][1])
        ob[i] = OS[row][col]
    return ob

#use OMP method to calculate initial coeffient values
def OMP(ob):
    rk = ob
    err = 0.0001*np.linalg.norm(ob,ord=2)
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
        mat[x][y] = 1
        u = fft.idct(fft.idct(mat,type=2,axis=-1,norm='ortho'),type=2,axis=-2,norm='ortho')
        for i in range(t):
            row = int(Position[i][0])
            col = int(Position[i][1])
            A[i][j] = u[row][col]
        B[:,j] = A[:,j]/np.linalg.norm(A[:,j])

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
        if ctn<err:
            xr = np.zeros(L)
            for i in range(k+1):
                xr[Skpos[i]] = xk[i]
            return xr

#calculate sampling points values
def getPointsValue(Coe1,Coe2,Coe3):
    point_val1 = np.zeros(t)
    point_val2 = np.zeros(t)
    point_val3 = np.zeros(t)
    for i in range(t):
        k = poslist[i]
        for j in range(L):
            point_val1[i] = point_val1[i] + Coe1[j]*f[j][k]
            point_val2[i] = point_val2[i] + Coe2[j]*f[j][k]
            point_val3[i] = point_val3[i] + Coe3[j]*f[j][k]
    return point_val1,point_val2,point_val3

#calculate Residual error
def getRes(Coe):
    Coe1,Coe2,Coe3 = divide(Coe)
    point_val1,point_val2,point_val3 = getPointsValue(Coe1,Coe2,Coe3)
    res = 0
    CoeMat1 = translate(Coe1)
    CoeMat2 = translate(Coe2)
    CoeMat3 = translate(Coe3)
    mapValue1 = fft.idct(fft.idct(CoeMat1,axis=-1,norm='ortho'),axis=-2,norm='ortho')
    mapValue2 = fft.idct(fft.idct(CoeMat2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
    mapValue3 = fft.idct(fft.idct(CoeMat3,axis=-1,norm='ortho'),axis=-2,norm='ortho')
    for i in range(t):
        res = res + (ob1[i]-point_val1[i])**2 + (ob2[i]-point_val2[i])**2 + (ob3[i]-point_val3[i])**2
    for i in range(1,M-1):
        for j in range(1,N-1):
            res = res + gamma*getPointRes(mapValue1,mapValue2[i][j],mapValue3,i,j)
    return res

def getPointRes(mat1,point_val,mat3,i,j):
    list1 = np.array([mat1[i-1][j-1],mat1[i-1][j],mat1[i-1][j+1],mat1[i][j-1],mat1[i][j],mat1[i][j+1],mat1[i+1][j-1],mat1[i+1][j],mat1[i+1][j+1]])
    list3 = np.array([mat3[i-1][j-1],mat3[i-1][j],mat3[i-1][j+1],mat3[i][j-1],mat3[i][j],mat3[i][j+1],mat3[i+1][j-1],mat3[i+1][j],mat3[i+1][j+1]])
    max1 = np.max(list1)
    min1 = np.min(list1)
    max3 = np.max(list3)
    min3 = np.min(list3)
    if max1 > point_val and min1 < point_val:
        point_res1 = 0
    if max1 <= point_val:
        point_res1 = (point_val-max1)**2
    if min1 >= point_val:
        point_res1 = (point_val-min1)**2
    if max3 > point_val and min3 < point_val:
        point_res3 = 0
    if max3 <= point_val:
        point_res3 = (point_val-max3)**2
    if min3 >= point_val:
        point_res3 = (point_val-min3)**2
    point_res = point_res1 + point_res3
    return point_res


#calculate gravity
def getGrav(Coe):
    Coe1,Coe2,Coe3 = divide(Coe)
    pRpc1 = np.zeros(L)
    pRpc2 = np.zeros(L)
    pRpc3 = np.zeros(L)
    pRpc = np.zeros(3*L)
    point_val1,point_val2,point_val3 = getPointsValue(Coe1,Coe2,Coe3)
        
    for j in range(L):
        for i in range(t):
            k = poslist[i]
            pRpc1[j] = pRpc1[j] + 2*f[j][k]*(point_val1[i]-ob1[i])
            pRpc2[j] = pRpc2[j] + 2*f[j][k]*(point_val2[i]-ob2[i])
            pRpc3[j] = pRpc3[j] + 2*f[j][k]*(point_val3[i]-ob3[i])
        
    for row in range(M):
        for col in range(N):
            
    
    return pRpc

#Barzilai-Borwein method
def BB(Coe1,Coe2,Coe3):
    iters = 0
    iter_max = 1000
    pre_m = 10
    eps = 10**-18
    beta = 0.5
    c1 = 0.1
    alpha = 0.0001
    alpha_max = 10000
    alpha_min = 0.01
    maxf = 0
    num = np.linspace(1,iter_max,iter_max)
    pre_res = np.zeros(iter_max)
    pre_alpha = np.zeros(iter_max)
    pre_pRpc = np.zeros(iter_max)
    tried_declines = np.zeros(iter_max)
    Coe = np.concatenate([Coe1,Coe2,Coe3])
    pRpc = getGrav(Coe)
    res = getRes(Coe)
    pre_pRpc[iters] = np.linalg.norm(pRpc,ord=2)
    plt.figure(figsize=(6, 6))
    while res > eps:
        tmp = c1*pre_pRpc[iters]**2
        pre_res[iters] = res
        count = min(iters,pre_m)
        tried_declines = pre_res - tmp*alpha
        res = getRes(Coe - alpha*pRpc)
        while res >= max(tried_declines[iters-count:iters+1]):
            alpha = beta*alpha
            tried_declines = pre_res - tmp*alpha
            new_Coe = Coe - alpha*pRpc
            res = getRes(new_Coe)
        pre_alpha[iters] = alpha
        Coe_pre = Coe
        Coe = Coe - alpha*pRpc
        pRpc_pre = pRpc
        pRpc = getGrav(Coe)
        sk = -alpha*pRpc_pre
        yk = pRpc - pRpc_pre
        s1 = np.dot(sk.T,yk)
        s2 = np.dot(yk.T,yk)
        alpha = s1/s2 #BB1
        #alpha = np.dot(sk.T,sk)/np.dot(sk.T,yk) #BB2
        if alpha > alpha_max:
            alpha = alpha_max
        if alpha < alpha_min:
            alpha = alpha_min
        iters = iters + 1
        pre_pRpc[iters] = np.linalg.norm(pRpc,ord=2)
        plt.semilogy(num[0:iters],pre_res[0:iters])
        plt.pause(0.01)
        if iters == 300:
            break
    Coe1,Coe2,Coe3 = divide(Coe)
    return Coe1,Coe2,Coe3

#calculate basis function
def getBasicFunction():
    f = np.zeros((L,L))
    for j in range(L):
        x=j//N
        y=j%N
        mat = np.zeros((M,N))
        mat[x][y] = 1
        u = fft.idct(fft.idct(mat,type=2,axis=-1,norm='ortho'),type=2,axis=-2,norm='ortho')
        for i in range(L):
            x2 = i//N
            y2 = i%N
            f[j][i] = u[x2][y2]
    return f
    
#translate coeffients to matrix format
def list2mat(Coe):
    CoeMat = np.zeros((M,N))
    for i in range(M):
        for j in range(N):
            idx = i*N + j
            CoeMat[i][j] = Coe[idx]
    return CoeMat

def mat2list(mat):
    Coe = np.zeros(M*N)
    for i in range(M):
        for j in range(N):
            idx = i*N + j
            Coe[idx] = mat[i][j]
    return Coe

def divide(Coe):
    return Coe[0:L],Coe[L:2*L],Coe[2*L:3*L]

#load parameters
inputfile1 = sys.argv[1]
inputfile2 = sys.argv[2]
inputfile3 = sys.argv[3]
outputfile1 = sys.argv[4]
outputfile2 = sys.argv[5]
outputfile3 = sys.argv[6]
pos = sys.argv[7]
OS1 = np.load(inputfile1,allow_pickle=True)
OS2 = np.load(inputfile2,allow_pickle=True)
OS3 = np.load(inputfile3,allow_pickle=True)
Position = np.load(pos,allow_pickle=True)

xfremax = 23
yfremax = 30
M = np.size(OS1,0)
N = np.size(OS1,1)   
t = Position.shape[0]
L = M*N
gamma = 0.00001
mu = 0.00001

poslist = np.zeros(t)
for i in range(t):
    poslist[i] = Position[i][0]*N + Position[i][1]

ob1 = getOriginalSamplingPointsValue(OS1)
ob2 = getOriginalSamplingPointsValue(OS2)
ob3 = getOriginalSamplingPointsValue(OS3)
f = getBasicFunction()
Coe1 = OMP(ob1)
Coe2 = OMP(ob2)
Coe3 = OMP(ob3)
StCoe1 = mat2list(fft.dct(fft.dct(OS1,type=2,axis=-1,norm='ortho'),type=2,axis=-2,norm='ortho'))
StCoe2 = mat2list(fft.dct(fft.dct(OS2,type=2,axis=-1,norm='ortho'),type=2,axis=-2,norm='ortho'))
StCoe3 = mat2list(fft.dct(fft.dct(OS3,type=2,axis=-1,norm='ortho'),type=2,axis=-2,norm='ortho'))
StCoe = np.concatenate([StCoe1,StCoe2,StCoe3])
StRes = getRes(StCoe)
Coe1,Coe2,Coe3 = BB(Coe1,Coe2,Coe3)
CoeMat1 = translate(Coe1)
CoeMat2 = translate(Coe2)
CoeMat3 = translate(Coe3)
RS1 = fft.idct(fft.idct(CoeMat1,axis=-1,norm='ortho'),axis=-2,norm='ortho')
RS2 = fft.idct(fft.idct(CoeMat2,axis=-1,norm='ortho'),axis=-2,norm='ortho')
RS3 = fft.idct(fft.idct(CoeMat3,axis=-1,norm='ortho'),axis=-2,norm='ortho')
plt.imshow(RS1)
plt.savefig('RS1.jpg')
plt.imshow(RS2)
plt.savefig('RS2.jpg')
plt.imshow(RS3)
plt.savefig('RS3.jpg')
np.save(outputfile1,RS1)
np.save(outputfile2,RS2)
np.save(outputfile3,RS3)
