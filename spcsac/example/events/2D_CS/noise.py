import numpy as np
import random
import matplotlib.pyplot as plt

def Noise(Mat):
    M,N = Mat.shape
    a = 0.000001
    b = -0.0001
    c = 0.01
    d = 0.015
    e = 0.000001
    f = -0.0003
    g = 0.01
    for i in range(M):
        for j in range(N):
            coe = 1 + a*i**3 + b*i**2 + c*i + d + e*j**3 + f*j**2 + g*j
            Mat[i][j] = coe*Mat[i][j]
    return Mat
