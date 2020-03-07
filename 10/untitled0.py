# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 20:58:56 2016

@author: Admin
"""
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import linalg as lin
from scipy.optimize import fsolve
from scipy.linalg import solve
from scipy.linalg import solve_banded

omega = 0.2
#lamb = 10
k = omega**2
sigma = 1/20
k0 = 50*np.pi
lamb = 0.25

def anal(x,t):
    alfa = k**(0.25)
    ksi = alfa*x
    ksil = alfa*lamb
    A = (alfa/(np.pi**0.5))**0.5
    return A*np.exp(-0.5*(ksi-ksil*np.cos(omega*t))**2 - 1j*(omega*t/2+ksi*ksil*np.sin(omega*t)-0.25*(ksil**2)*np.sin(2*omega*t)))
    
def zac(x):
    alfa = k**(0.25)
    A = (alfa/(np.pi**0.5))**0.5
    return A*np.exp(-alfa**2 *0.5*(x-lamb)**2)

def zac2(x):
    A = (2*np.pi*sigma**2)**(-0.25)
    return A*np.exp(1j*k0*(x-lamb))*np.exp(-(x-lamb)**2/(4*sigma**2))    
def anal2(x,t):
    A = (2*np.pi*sigma**2)**(-0.25)/((1+1j*t/(2*sigma**2))**0.5)
    return A*np.exp((-(x-lamb)**2/(4*sigma**2)+1j*k0*(x-lamb) - 1j*k0**2*t/2)/(1+1j*t/(2*sigma**2)))

def naredimA(N,a,b,V,dt):
    mA = []
    for i in range(1,N-1):
        vrstica = []
        for j in range(1,N-1):    
            if i==1:
                if j==1:
                    vrstica.append(d(1,b,V,dt))
                elif j==2:
                    vrstica.append(a)
                else:
                    vrstica.append(0)
            elif i==(N-2):
                if j==(N-2):
                    vrstica.append(d(N-2,b,V,dt))
                elif j==(N-3):
                    vrstica.append(a)
                else:
                    vrstica.append(0)
            else:
                if i==j:
                    vrstica.append(d(i,b,V,dt))
                elif abs(i-j) == 1:
                    vrstica.append(a)
                else:
                    vrstica.append(0)
        mA.append(vrstica)
    return np.matrix(mA)

def naredimAband(N,a,b,V,dt):
    mAband = []
    aji = [a for i in range(1,N-2)]
    prva = aji[:]
    prva.insert(0,0)
    tretja = aji[:]
    tretja.append(0)
    druga = [d(j,b,V,dt) for j in range(1,N-1)]
    mAband.append(prva)
    mAband.append(druga)
    mAband.append(tretja)
    return mAband    


def d(j,b,V,dt):
    #return 1+b+1j*0.5*dt*V[j]
    return 1+b

def zracunaj(cas,stevilo,stevilot):
    x = np.linspace(-0.5,1.5,stevilo)
    t = np.linspace(0,cas,stevilot)
    N = len(x)
    dx = (x[-1]-x[0])/N
    NN = len(t)
    dt = (t[-1]-t[0])/NN
    psi0= zac2(x)
    psi0 = psi0.tolist()
    psi0=psi0[1:-1] #prvo in zadnje je 0
    psi0 = np.matrix(psi0).T
    b = 1j * 0.5*dt/(dx**2)
    a = -b/2
    V = 0.5*k*(x**2)
    #V = [0.5*k*(x[i]**2) for i in range(len(x))]
    mA = naredimA(N,a,b,V,dt)
    mAband = naredimAband(N,a,b,V,dt)
    mA = mA.H.T
    for i in t:
        desna = mA*psi0
        #print(mAband)
        #print(desna)
        psi0 = solve_banded((1,1),mAband,desna)
    psi0 = psi0.tolist()
    psi0.insert(0,[0])
    psi0.append([0])
    psi0 = np.array(psi0)
    return [x,psi0]
    
#V = V[1:-1] #na robu je nič
    

def napaka(cas,N,NN):
    x = np.linspace(-0.5,1.5,N)
    t = np.linspace(0,cas,NN)
    prava = anal2(x,cas)
    prava = abs(prava)**2
    aproks = zracunaj(cas,N,NN)[1]
    aproks = abs(aproks)**2
    aproks = aproks.flatten()
    napaka = abs(prava-aproks)
    return [x,napaka]

def napaka2(cas,N,NN):
    x = np.linspace(-40,40,N)
    t = np.linspace(0,cas,NN)
    prava = anal(x,cas)
    prava = abs(prava)**2
    aproks = zracunaj(cas,N,NN)[1]
    aproks = abs(aproks)**2
    aproks = aproks.flatten()
    napaka = abs(prava-aproks)/abs(prava)
    return [x,napaka]

per = 2*np.pi /omega
"""
x = np.linspace(-40,40,500)
y = abs(anal(x,per))**2
plt.plot(x,y,label="t=T")
y = abs(anal(x,5*per/4))**2
plt.plot(x,y,label="t=5T/4")
y = abs(anal(x,3*per/2))**2
plt.plot(x,y,label="t=3T/2")
y = abs(anal(x,7*per/4))**2
plt.plot(x,y,label="t=7T/4")
y = abs(anal(x,2*per))**2
plt.plot(x,y,label="t=2T")
plt.legend()
plt.title("Analitična rešitev")
plt.savefig("nihaji2a.pdf")
plt.show()
"""



#rez = zracunaj(3*per,300,1000)
#plt.plot(rez[0],abs(rez[1])**2,label="t=3T")
#rez = zracunaj(4*per,300,1000)
#plt.plot(rez[0],abs(rez[1])**2,label="t=4T")
#rez = zracunaj(5*per,300,1000)
#plt.plot(rez[0],abs(rez[1])**2,label="t=5T")
#rez = zracunaj(8*per,300,1000)
#plt.plot(rez[0],abs(rez[1])**2,label="t=8T")
"""
sa = 0.0005
deltat = sa/1000
deltax = (deltat/2)**0.5
Nx = int(2/deltax)
rez = zracunaj(sa,Nx,1000)
plt.plot(rez[0],abs(rez[1])**2,label=r"$t = 0.0005$")
"""
"""
x = np.linspace(-0.5,1.5,100)
plt.plot(x,abs(zac2(x))**2,label=r"$t=0$")
sa = 0.0005
deltat = sa/1000
deltax = (deltat/2)**0.5
Nx = int(2/deltax)
rez = zracunaj(sa,Nx,1000)
plt.plot(rez[0],abs(rez[1])**2,label=r"$t = 0.0005$")
sa = 0.001
deltat = sa/1000
deltax = (deltat/2)**0.5
Nx = int(2/deltax)
rez = zracunaj(sa,Nx,1000)
plt.plot(rez[0],abs(rez[1])**2,label=r"$t = 0.001$")
sa = 0.0015
deltat = sa/1000
deltax = (deltat/2)**0.5
Nx = int(2/deltax)
rez = zracunaj(sa,Nx,1000)
plt.plot(rez[0],abs(rez[1])**2,label=r"$t = 0.0015$")
sa = 0.002
deltat = sa/1000
deltax = (deltat/2)**0.5
Nx = int(2/deltax)
rez = zracunaj(sa,Nx,1000)
plt.plot(rez[0],abs(rez[1])**2,label=r"$t = 0.002$")
sa = 0.0025
deltat = sa/1000
deltax = (deltat/2)**0.5
Nx = int(2/deltax)
rez = zracunaj(sa,Nx,1000)
plt.plot(rez[0],abs(rez[1])**2,label=r"$t = 0.0025$")
sa = 0.003
deltat = sa/1000
deltax = (deltat/2)**0.5
Nx = int(2/deltax)
rez = zracunaj(sa,Nx,1000)
plt.plot(rez[0],abs(rez[1])**2,label=r"$t = 0.003$")
sa = 0.0035
deltat = sa/1000
deltax = (deltat/2)**0.5
Nx = int(2/deltax)
rez = zracunaj(sa,Nx,1000)
plt.plot(rez[0],abs(rez[1])**2,label=r"$t = 0.0035$")
plt.legend()
plt.title(r"$N_t=1000$")
plt.savefig("Druga/potek3.pdf")
plt.show()
"""
"""
x = np.linspace(-40,40,300)
y = anal(x,100)
plt.plot(x,abs(y)**2)
a = zracunaj(37,300,1000)
plt.plot(a[0],abs(a[1])**2)
plt.show()
"""

"""
n = napaka(0.0001,300,1000)
plt.plot(n[0],n[1],label="")
plt.legend()
plt.title(r" Absolutna napaka aproksimacije za več $N_t,N_x$")
plt.savefig("Druga/1.pdf")
plt.show()
"""

"""
print(per)
x = np.linspace(-40,40,100)
y = np.linspace(0,per,10)
z=[]
for i in y:
    ars = zracunaj(i,300,1000)[1]
    #print(np.array(ars).flatten())
    ars = np.array(ars).flatten()
    z.append(list(ars))
x, y = np.meshgrid(x,y)
fig = plt.figure()
ax = fig.gca(projection="3d")
plt.xlabel("x")
plt.ylabel("t")
plt.title(r"Aproksimacija, 1 nihaj")
ax.plot_surface(x,y,z,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0)
plt.savefig("aproks.pdf")
plt.show()
"""

x = np.linspace(-0.5,1.5,1000)
y = abs(anal2(x,0.000))**2
plt.plot(x,y,label="t=0")
y = abs(anal2(x,0.0005))**2
plt.plot(x,y,label="t=0.0005")
y = abs(anal2(x,0.001))**2
plt.plot(x,y,label="t=0.001")
y = abs(anal2(x,0.0015))**2
plt.plot(x,y,label="t=0.0015")
y = abs(anal2(x,0.002))**2
plt.plot(x,y,label="t=0.002")
y = abs(anal2(x,0.0025))**2
plt.plot(x,y,label="t=0.0025")
y = abs(anal2(x,0.003))**2
plt.plot(x,y,label="t=0.003")
y = abs(anal2(x,0.0035))**2
plt.plot(x,y,label="t=0.0035")
plt.legend()
plt.title("Analitična rešitev")
plt.savefig("Druga/anal.pdf")
plt.show()    
    
    
    
    
    
    
    
    
def PsiT_2(x, t, s=0.05, l=0.25, k0=50*np.pi):
    koef1 = 1/np.sqrt(np.sqrt(2*np.pi*s**2))
    koef2 = 1/np.sqrt(1 + 0.5*1j*t/s**2)
    zgoraj = -0.25*(x - l)**2/s**2 + 1j*k0*(x-l) - 0.5*1j*t*k0**2
    spodaj = 1 + 0.5*1j*t/s**2
    return koef1*koef2 * np.exp(zgoraj/spodaj)    
    
    
    
    
    
    
    
    
    
    
    