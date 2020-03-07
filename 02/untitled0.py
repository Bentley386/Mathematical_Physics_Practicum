# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 00:16:36 2016

@author: Admin
"""

import matplotlib.pyplot as plt
import numpy as np
from random import *
from scipy.optimize import curve_fit


def kot():
    m = 2*np.pi
    return random()*m
def dolzina(mu):
    r = random()
    return (1-r)**(1/(1-mu))


def koraki(n,mu,stop=0):    
    x = [0]
    y = [0]
    lji = []
    cas = 0
    for i in range(n):
        l = dolzina(mu)
        cas+=l
        a = kot()
        x.append(x[-1]+l*np.cos(a))
        y.append(y[-1]+l*np.sin(a))
        lji.append(cas)
        if stop!=0 and cas>=stop:
            break
    return [x,y,lji]

def koraki2(n,mu,stop=0):    
    x = [0]
    y = [0]
    lji = []
    cas = 0
    for i in range(n):
        xi = abs(random()-0.5)
        x.append(x[-1]+xi)
        yi = abs(random()-0.5)
        y.append(y[-1]+yi)
        lji.append((x[-1]**2+y[-1]**2)**0.5)
        cas+=((xi**2+yi**2)**0.5)
        if stop!=0 and cas>=stop:
            break
    return [x,y,lji]

def koraki3(n,mu,stop=0,nu=1.5):    
    x = [0]
    y = [0]
    lji = []
    cas = 0
    for i in range(n):
        cas+=dolzina(nu)
        l = dolzina(mu)
        a = kot()
        x.append(x[-1]+l*np.cos(a))
        y.append(y[-1]+l*np.sin(a))
        cas+=l
        lji.append(cas)
        if stop!=0 and cas>=stop:
            break
    return [x,y,lji]
def mad(x):
    x = np.array(x)
    return np.median(abs(x-np.median(x)))

def deviacija(t,mu):
    """izracuna deviacijo pri nekem casu"""
    rji = []
    temp = []
    for i in range(800):
        temp = koraki3(1000000,mu,t)
        rji.append((temp[0][-1]**2+temp[1][-1]**2)**0.5)
    return mad(rji)

def f(x,a):
    """nastavek za varianco"""
    return x**(a/2)

prvo = koraki(100000,2.2,10000)
prvok = range(1,len(prvo[2])+1)
drugo = koraki3(10000,2.2,10000)
drugok=range(1,len(drugo[2])+1)
    

    
#parametri = curve_fit(f,[100,500,1000,5000,10000],[4.2,9.2,11.8,26.1,33.6])
#print(parametri[0])
"""
#t = np.linspace(0,100,1000)
plt.plot(prvo[2],prvok,"k",label=r"Brez sticking time") #x
plt.plot(drugo[2],drugok,"b",label=r"Z sticking time ($\nu = 1.5$)") #x abs
#plt.plot(t,t**1.6,"g",label=r"$f(x)=\frac{1}{x}$") #1/x
#plt.plot(t,t**1.25,"r",label=r"$f(x)=\frac{1}{|x|}$")  #1/absx
#plt.plot(t,t**2.1,"c",label=r"$\mu = 1.5$") #mu 1.5
#plt.plot(t,t**1.16,"m",label=r"$\mu=2.5$") #mu 2.5
#plt.plot(t,t**1.4,"y",label=r"$\mu=2.2$") #mu 2.2
plt.legend(loc="upper left")
plt.title(r"Število sprehodov narejenih po določenem času, $\mu=2.2$")
#plt.yscale("log")
#plt.xscale("log")
plt.savefig("zadnjeres.pdf")
plt.show()
"""
print(10000**0.38)
"""
rezultat = koraki(10000,2.2,10000)
print(len(rezultat[0]))
plt.plot(rezultat[0],rezultat[1],"k.")
plt.plot(rezultat[0],rezultat[1],"k")
plt.title(r"$\mu = 2.2,\nu = 1.5 , t=10000$ Z sticking time")
plt.savefig("z.pdf")
plt.show()
"""
"""
l = koraki(10000,2.5)[2]
m = koraki(10000,2.2)[2]
n = koraki(10000,1.5)[2]

plt.hist([l,m,n],bins=[1,2,3,4,5,6,7,8,9])
plt.title(r"$Modra:\mu=2.5, Zelena:\mu=2.2, Rdeča: \mu=1.5$")
plt.savefig("histogram.pdf")
"""