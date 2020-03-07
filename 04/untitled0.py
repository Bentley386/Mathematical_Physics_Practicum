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

def disk(funk,N,interval,hann=False):#naredi diskretne vrednosti funkcije
    x = []
    delta = (interval[1]-interval[0])/N
    vrednosti = []
    for i in range(N):
        x.append(i*delta)
        aaa  = funk(i*delta)
        if hann==True:
            aaa = aaa*0.5*(1-np.cos(2*np.pi*i/(N-1))) #window function
        vrednosti.append(aaa)
    return [x,vrednosti]
    
def Hn(vred,n): #zracuna Hn od nekih vrednosti
    N = len(vred)
    suma=0
    for i in range(N):
        suma+= vred[i]*np.exp(2*np.pi*1j*i*n/N)
    return suma
def hk(vred,n):
    N = len(vred)
    suma=0
    for i in range(N):
        suma+= vred[i]*np.exp(-2*np.pi*1j*i*n/N)
    return suma/N
    
def sortiraj(fourier,N,interval):
    koef = []
    delta = (interval[1]-interval[0])/N
    for i in range(len(fourier)):
        koef.append(fourier[i]/delta)
    koef = koef[:-1] #vzamemo proc zadnjega ker se itak ponovi na prvem
    sredina = N//2
    x = []
    y = []
    i = 0
    j = sredina
    k = 0
    while(len(x)!=N):
        x.append(delta*i)
        if j<len(koef):
            y.append(koef[j])
        else:
            y.append(koef[k])
            k+=1
        j+=1
        i+=1
    return [x,y]
        


def inverz(fourier,N,interval):
    y=[]
    vred = sortiraj(fourier[1],N,interval)
    for i in range(N):
        y.append(hk(vred[1],i))  
    return [vred[0],y]
        
        
def fourier(vred,N,interval): #vrnex in y v fourierevem prostoru od vrednosti
    x = []
    y = []
    delta = (interval[1]-interval[0])/N
    n=-N/2
    while(n<=N/2):
        x.append(n/(N*delta))
        y.append(delta*Hn(vred,n))
        n+=1
    return [x,y]
    
    
def moc(four): #iz fourier transf naredi moc
    x=[]
    moci=[]
    N=len(four[0])
    for i in range(N):
        if four[0][i]<0:
            continue
        elif four[0][i]==0:
            x.append(four[0][i])
            moci.append(abs(four[1][i])**2)
        elif four[0][i]>0:
            x.append(four[0][i])
            moci.append(2*abs(four[1][i])**2)
    return [x,moci]


def prva(x):
    return 3*np.sin(3*np.pi*x)+2*np.sin(2.5*np.pi*x)
def druga(x):
    return np.sin(2*np.pi*x)
    
def okni(x):
    y = []
    N = len(x)
    for i in range(N):
        aaa = x[i]*0.5*(1-np.cos(2*np.pi*i/(N-1))) #window function
        y.append(aaa)
    return y

"""with open("akus/Previden.txt","r") as fajl:
    teks=fajl.read()
teks = teks.split()
teks = list(map(int,teks))
"""

with open("akus/zelomocen.dat","r") as fajl: #frekv 44.1kHz, 0-1s int
    teks = fajl.read()
    teks = teks.split()
    teks = list(map(float,teks))
        
print(teks)
print(len(teks))  
    
#vrednosti = disk(prva,5000,[0,8.01])
f = fourier(okni(teks),len(teks),[0,1])

print("ja")
#interval je 0 do 1.325
#vrednosti= disk(prva,1000,[0,5])
#plt.plot(vrednosti[0],vrednosti[1])
#plt.xlabel("t")
#plt.ylabel("h(t)")
#plt.title(r"$h(t) = 3sin(3\pi t) + 2sin(2.5\pi t)$, N=1000")
#plt.savefig("druga/sin.pdf")
#plt.show()



#inv = inverz(f,1000,[0,8.01])
#plt.plot(inv[0],inv[1],"b",label="Inverz transformacije")
#plt.plot(vrednosti[0],vrednosti[1],"r",label="Original signal")
#plt.legend()
#plt.xlabel("t")
#plt.ylabel("h(t)")
#plt.savefig("prva/inverz.pdf")
#plt.show()
moci = moc(f)
plt.plot(moci[0],moci[1],"b",label="Brez oknenja")
plt.yscale("log")
plt.title(r"$P(\nu)$, Akus. reso. zelo močen udarec, 44.1kHz")
plt.xlabel(r"$\nu$")
plt.ylabel(r"$P(\nu)$")
#vrednosti = disk(prva,5000,[0,8.01],True)
#f = fourier(vrednosti[1],5000,[0,8.01])
#moci=moc(f)
#plt.plot(moci[0],moci[1],"r",label="Z oknenjem")
#plt.legend(loc="lower left")
#plt.xlim(0,4)
#plt.ylim(0.01,13)

#plt.plot(moci[0],moci[1],"bo",label="Brez oknenja")
#vrednosti = disk(druga,1000,[0,4.66],True)
#f = fourier(vrednosti[1],1000,[0,4.66])
#moci = moc(f)
#plt.plot(moci[0],moci[1],"ro",label="Z oknjenjem")
#plt.legend(loc="lower right")
#plt.title("h(t) N=1000")
#plt.xlim(0,2)
plt.savefig("akus/zelo1.pdf")
plt.clf()
plt.plot(moci[0],moci[1],"b",label="Brez oknenja")
plt.yscale("log")
plt.title(r"$P(\nu)$, Akus. reso. zelo močen udarec, 44.1kHz")
plt.xlabel(r"$\nu$")
plt.ylabel(r"$P(\nu)$")
plt.xlim(0,500)
plt.savefig("akus/zelo2.pdf")
plt.clf()
plt.plot(moci[0],moci[1],"b",label="Brez oknenja")
plt.yscale("log")
plt.title(r"$P(\nu)$, Akus. reso. zelo močen udarec, 44.1kHz")
plt.xlabel(r"$\nu$")
plt.ylabel(r"$P(\nu)$")
plt.xlim(0,1000)

plt.savefig("akus/zelo3.pdf")
plt.clf()
plt.plot(moci[0],moci[1],"b",label="Brez oknenja")
plt.yscale("log")
plt.title(r"$P(\nu)$, Akus. reso. zelo močen udarec, 44.1kHz")
plt.xlabel(r"$\nu$")
plt.ylabel(r"$P(\nu)$")
plt.xlim(500,1000)
plt.savefig("akus/zelo4.pdf")
plt.show()
