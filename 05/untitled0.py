# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:36:40 2016

@author: Admin
"""

from scipy.fftpack import fft
import numpy as np
from scipy.fftpack import ifft
import matplotlib.pyplot as plt

def okni(x):
    y = []
    N = len(x)
    for i in range(N):
        aaa = x[i]*0.5*(1-np.cos(2*np.pi*i/(N-1))) #window function
        y.append(aaa)
    return y



def kvadriraj(x):
    for i in range(len(x)):
        x[i] = abs(x[i])**2
    return x

def zmnozi(x,y):
    z = []
    for i in range(len(x)):
        z.append(x[i]*np.conj(y[i]))
    return z
def sortiraj(z):
    #z je input
    N = len(z)
    n = -N/2+1 #x koordinate
    x = []
    y = []
    i = N//2  #sredina
    j = 1
    while len(x)!=(N-1):
        x.append(n)
        n+=1
        if i<N:
            y.append(z[i])
            i+=1
        else:
            y.append(z[j])
            j+=1
    return [x,y]

def reskaliraj(z):
    h = (sum(a)/2**11)**2
    nic = z[len(z)//2]
    aa = nic-h
    for i in range(len(z)):
        z[i] = (z[i]-h)/aa
    return z
    
def reskaliraj2(z):
    h=sum(map(lambda x: x**2,a))
    n = len(a)
    for i in range(len(z)):
        z[i] = z[i]/h
    return z
def reskaliraj3(z):
    h = abs(sum([a[i]*b[i] for i in range(len(a))]))
    for i in range(len(z)):
        z[i] = z[i]/h
    return z
def reskaliraj4(z):
    for i in range(len(a)):
        z[i] = z[i]/abs(sum(a)*sum(b))
    return z
    
def reskaliraj5(z):  #ta je prava...
    f = sum(list(map(lambda x: x**2,a)))**0.5
    g = sum(list(map(lambda x : x**2,a)))**0.5
    for i in range(len(z)):
        z[i] = z[i]/(f*g)
    return z
    
def reskaliraj6(z):
    f = sum(list(map(lambda x: x**2,a)))
    for i in range(len(z)):
        z[i] = z[i]/(f)
    return z
with open("10m_radioactivitymix.txt","r") as f:
    a = f.read()
    a = list(map(float,a.split()))
    print(len(a))
"""
with open("cricki.txt","r") as f:
    b = f.read()
    b = list(map(float,b.split()))
    print(len(b))
"""
def normaliziraj(z):
    for i in range(len(a)):
        z[i] = z[i]/(abs(a[i]))**2
    return z
#auta = okni(a) oknenje pokvari

def dopolni(x,a):
    while len(x)!=a:
        x.append(0)
    return x
print(2**13)
fourier1 = fft(a,2**13)
#fourier2 = fft(b,2**19)

korelac = sortiraj(ifft(kvadriraj(fourier1)))
y = reskaliraj6(np.real(korelac[1]))
#y = np.real(korelac[1])
plt.plot(korelac[0],y)
#plt.xlim(-1000,1000)
plt.title("Radioaktivnost zemlje s periodično komponento")
plt.ylabel("Avtokorelacija")
plt.xlabel("Zamik v enotah 10min")
plt.savefig("radioo.pdf")
plt.show()

"""
x = np.linspace(1700,len(a),len(a))
plt.plot(x,a,"b.")
plt.title("Radioaktivnost zemlje s periodično komponento")
plt.savefig("radioo2.pdf")
plt.show()
"""