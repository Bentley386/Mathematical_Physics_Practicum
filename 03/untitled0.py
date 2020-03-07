# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 01:39:12 2016

@author: Admin
"""
import matplotlib.pyplot as plt
from scipy.misc import factorial
from scipy.special import hermite
import numpy as np
from random import *
from numpy.linalg import eigh #da not matriko z arrayi, vrne prvo eigenvalues ki narascajo in
#druog matriko z eigenvectors

def disk(funk,delta,N):
    vrednosti = []
    for i in range(N):
        vrednosti.append(funk(i*delta))
    return vrednosti
    
def fourier(vred):
    N = len(vred)
    suma=0
    for i in range(N):
        suma+= vred[i]*np.exp(2*np.pi*j*i*n/N)

print(1+5j+2j)