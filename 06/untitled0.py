# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 15:28:20 2016

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.optimize import fmin
from scipy.optimize import fsolve
from scipy.optimize import bisect
def prava(t0,t,h,k):
    """Vrne x in y x je čas y pa vrednosti od t=0 do časa t"""
    x = [0]
    y = [t0]
    cas = 0
    while cas<t:
        cas+=h
        x.append(cas)
        ips = -5 + np.exp(-k*cas)*(t0 + 5)
        y.append(ips)
    return [x,np.array(y)]

def euler(t0,t,h,k):
    x = [0]
    y = [t0]
    cas = 0
    while cas<t:
        cas+=h
        x.append(cas)
        ips = y[-1] - h*k*(y[-1]+5)
        y.append(ips)
    return [x,np.array(y)]

            
def euler2(t0,t,h,k):
    x = [0,h]
    y = [t0,-5+np.exp(-k*h)*(t0+5)]
    cas = h
    while cas<t:
        cas+=h
        x.append(cas)
        ips = -2*h*k*(y[-1]+5)+y[-2]
        y.append(ips)
    return [x,np.array(y)]
            
def euler3(t0,t,h,k):
    x = [0]
    y = [t0]
    cas = 0
    while cas<t:
        cas+=h
        x.append(cas)
        AA = 1+k*h/2
        ips = (y[-1]*(1-k*h/2) - k*h*5)/AA
        y.append(ips)
    return [x,np.array(y)]

def rkp(x,y,A,z): #rk pomozna
    return -0.1*(y-z) + A*np.sin(2*np.pi/24*(x-10))
def euler4(t0,t,h,A,z): #da mi ni treba vse spremenit v rk4
    x = [0]
    y = [t0]
    cas = 0
    while cas<t:
        cas+=h
        x.append(cas)
        k1 = rkp(x[-1],y[-1],A,z)
        k2 = rkp(x[-1]+0.5*h,y[-1]+0.5*h*k1,A,z)
        k3 = rkp(x[-1]+0.5*h,y[-1] + 0.5*h*k2,A,z)
        k4 = rkp(x[-1]+h,y[-1] + h*k3,A,z)
        ips = y[-1] + h/6*(k1+2*k2+2*k3+k4)
        y.append(ips)
    return [x,np.array(y)]

def okni(x):
    l = len(x)-1
    rezultat = []
    for i in range(len(x)):
        rezultat.append(x[i]*(np.sin(np.pi*i/l)**2))
    return np.array(rezultat)
[x,y] = euler4(21,100,0.001,1,-5)


"""yy = okni(y)            
fourier = fft(yy)            
def spekter(x,y):
    x1=[0]
    y1=[abs(y[0])**2]
    for i in range(1,len(x)//2):
        x1.append(i)
        y1.append(2*abs(y[i])**2)
    return [x1,y1]
[x1,y1] = spekter(x,fourier)               """

# = euler4(21,100,0.001,100)


#for i in range(len(x)):
    #y[i] = -0.1*(y[i]+5) + np.sin(2*np.pi/24*(x[i]-10))

def f(g):
    najblizja = (np.abs(np.array(x)-g)).argmin()
    return y[najblizja]    
#nicle = bisect(f,90,95)
#print(nicle)    
#========za vec h============
plt.plot(x,y)
plt.arrow(19.13,-5,0,100)
plt.arrow(44.5,-5,0,100)
plt.arrow(68.596,-5,0,100)
plt.arrow(92.6,-5,0,100)
#plt.plot(aa,np.sin(aa))
#plt.yscale("log")
plt.xlim(0,100) #max pri 42
#plt.plot(euler4(21,100,0.001,10,40)[0],euler4(21,100,0.001,10,40)[1],label=r"A=10")
#plt.plot(euler4(21,100,0.001,0.1,40)[0],euler4(21,100,0.001,0.1,40)[1],label=r"A=0.1")
#plt.plot(euler4(21,100,0.001,-5,40)[0],euler4(21,100,0.001,-5,40)[1],label=r"A=-5")



#plt.plot(euler4(21,100,0.1,10,-5)[0],euler4(21,100,0.1,10,-5)[1],label="h=0.1")
#plt.plot(euler4(21,100,5,10,-5)[0],euler4(21,100,5,10,-5)[1],label="h=5")
#plt.plot(euler4(21,100,0.01,10,-5)[0],euler4(21,100,0.01,10,-5)[1],label="h=0.01")
#========= relativna napaka==========
#plt.plot(prava(21,10,0.1,10)[0],abs(prava(21,10,0.1,10)[1]-euler4(21,10,0.1,10)[1])/abs(prava(21,10,0.1,10)[1]),label="h=0.1")
#plt.plot(euler4(21,100,0.001,1)[0],abs(euler4(21,100,0.001,1)[1]-euler4(21,100,0.1,1)[1])/abs(euler4(21,100,0.001,1)[1]),label="h=0.1,A=1")
#plt.plot(euler4(21,100,0.001,1)[0],abs(euler4(21,100,0.001,1)[1]-euler4(21,100,1,1)[1])/abs(euler4(21,100,0.001,1)[1]),label="h=1, A=1")
#plt.plot(euler4(21,100,0.001,1)[0],abs(euler4(21,100,0.001,1)[1]-euler4(21,100,5,1)[1])/abs(euler4(21,100,0.001,1)[1]),label="h=5, A=1")
#plt.plot(euler4(21,100,0.001,0)[0],abs(euler4(21,100,0.001,0)[1]-euler4(21,100,0.1,0)[1])/abs(euler4(21,100,0.001,0)[1]),label="h=0.1 ,A=0")
#plt.plot(euler4(21,100,0.001,0)[0],abs(euler4(21,100,0.001,0)[1]-euler4(21,100,1,0)[1])/abs(euler4(21,100,0.001,0)[1]),label="h=1, A=0")
#plt.plot(euler4(21,100,0.001,0)[0],abs(euler4(21,100,0.001,0)[1]-euler4(21,100,5,0)[1])/abs(euler4(21,100,0.001,0)[1]),label="h=5, A=0")

#========absolutna napaka==============
#plt.plot(euler4(21,10,0.001,100)[0],abs(euler4(21,10,0.001,100)[1]-prava(21,10,0.001,100)[1]),label="h=0.001")
#plt.plot(euler4(21,10,0.01,100)[0],abs(euler4(21,10,0.01,100)[1]-prava(21,10,0.01,100)[1]),label="h=0.01")
#plt.plot(euler4(21,10,0.005,100)[0],abs(euler4(21,10,0.005,100)[1]-prava(21,10,0.005,100)[1]),label="h=0.005")
#plt.plot(euler4(21,10,0.005,10)[0],abs(euler4(21,10,0.005,10)[1]-prava(21,10,0.005,10)[1]),label="h=0.005")
#plt.plot(euler4(21,400,5)[0],abs(euler4(21,400,5)[1]-prava(21,400,5)[1]),label="h=5")

#plt.legend(loc="lower right")
plt.xlabel("t")
plt.ylabel("T")
#plt.yscale("log")
plt.title(r"Maksimumi")
plt.savefig("dod/max3.pdf")
plt.show()            
            
            
            
            
            
            
            
            
            
            
            
"""
PRVA EULERJEVA
a = prava(21,100,0.1)
plt.plot(a[0],a[1],label="Prava odvisnost")
plt.plot(euler(21,100,20)[0],euler(21,100,20)[1],label="h=20")
#plt.plot(prava(21,100,0.1)[0],abs(prava(21,100,0.1)[1]-euler(21,100,0.1)[1])/abs(prava(21,100,0.1)[1]),label="h=0.1")
#plt.plot(prava(21,100,0.2)[0],abs(prava(21,100,0.2)[1]-euler(21,100,0.2)[1])/abs(prava(21,100,0.2)[1]),label="h=0.2")
#plt.plot(prava(21,100,0.5)[0],abs((prava(21,100,0.5)[1]-euler(21,100,0.5)[1])/(prava(21,100,0.5)[1])),label="h=0.5")
#plt.plot(prava(21,100,1)[0],abs(prava(21,100,1)[1]-euler(21,100,1)[1])/abs(prava(21,100,1)[1]),label="h=1")
plt.legend(loc="upper left")
plt.xlabel("t")
plt.ylabel("T")
plt.xlim(0,100)
#plt.yscale("log")
#plt.ylim(6,12)
plt.title("T(t),relativna razlika Eulerjevih metod od prave vrednosti za več h")
plt.savefig("euler/12.pdf")
plt.show()
"""