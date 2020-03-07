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
from scipy.special import jn_zeros #prvi parameter je order, drugi št. ničel
from scipy.special import jv #prvi order drugi argument


def zracun(m,s,k):
    C = 0
    for i in range(m):
        nicle = jn_zeros(2*i+1,s+5)
        for j in range(1,s):
            I = 0
            for l in range(i+1,k+i+1):
                I+= l*jv(2*l,nicle[j-1])/(4*l**2 -1)
            I = I*4*(2*i+1)/nicle[j-1]
            C+= (8/np.pi * I/(2*i+1)/(nicle[j-1]*jv(2*i+2,nicle[j-1])))**2
    C = 8*C
    return C
def vecm():
    x = [a for a in range(1,101)]
    y = []
    for i in x:
        y.append(zracun(10,20,i))
    return [x,y]
def vrednost(x,fi):
    suma =0
    for i in range(10):
        nicle = jn_zeros(2*i+1,45)
        for j in range(1,21):
            I = 0
            for l in range(i+1,i+11):
                I+= l*jv(2*l,nicle[j-1])/(4*l**2 -1)
            I = I*8/nicle[j-1] #modificiran I
            suma+= 4*I/(np.pi*jv(2*i+2,nicle[j-1])**2)*jv(2*i+1,nicle[j-1]*x)*np.sin((2*i+1)*fi)
    return suma
"""  
#m=10 s=10   
I = 0
i = 10
j = 20
nicle = jn_zeros(2*i+1,50)
for l in range(i+1,20+i+1):
    I+= l*jv(2*l,nicle[j-1])/(4*l**2 -1)
I = I*4*(2*i+1)/nicle[j-1] #modificiran I    
print(I)
"""    
    
def csuma(apb,n):
    suma = 0
    for i in range(1,n+1):
        suma+= np.tanh(apb*(2*i-1)*np.pi/2)/(2*i-1)**5
    return suma
def canal(apb,n):
    return abs(2*np.pi/apb * (1/3 - 64/(apb*np.pi**5))*csuma(apb,n))
def zracunaj(apb):
    x = [i for i in range(1,21)]
    y = [canal(apb,i) for i in x]
    return [x,y]
print(canal(1,100000))
"""
x = np.linspace(1,10)
y = canal(x,200)
plt.plot(x,y)
plt.title("Odvisnost C-ja od razmerja a/b")
plt.ylabel("C")
plt.xlabel("a/b")
plt.savefig("druga2.pdf")
plt.show()
"""
"""
a = zracunaj(1)
plt.plot(a[0],a[1],label=r"$\frac{a}{b}=1$")
a = zracunaj(2)
plt.plot(a[0],a[1],label=r"$\frac{a}{b}=2$")
a = zracunaj(0.5)
plt.plot(a[0],a[1],label=r"$\frac{a}{b}=0.5$")
a = zracunaj(5)
plt.plot(a[0],a[1],label=r"$\frac{a}{b}=5$")
a = zracunaj(10)
plt.plot(a[0],a[1],label=r"$\frac{a}{b}=10$")
plt.title("C za več razmerij a/b")
plt.xlabel("Število šeštetih členov")
plt.ylabel("C")
plt.legend()
plt.savefig("druga1.pdf")
plt.show()
"""

def lenaloga(m,n,apb):
    C = 0
    i = 1
    j = 1
    while i<=m:
        while j<=n:
            C+= ((i*j)**2*((j*np.pi)**2 /apb + apb*(i*np.pi)**2))**(-1)
            j+=2
        i+=2
    return 32*16/(np.pi**3) * C

def uiks(x,y,b):
    i = 1
    j= 1
    suma = 0
    while i<=11:
        while j<=11:
            koef = 16/(i*j*np.pi**2)/((i*np.pi)**2 +(j*np.pi/b)**2)
            suma += koef* np.sin(i*np.pi*x) * np.sin(j*np.pi * y/b)
            j+=2
        i+=2
    return suma
"""    
b = 0.2
x = np.linspace(0,1,200)    
y = np.linspace(0,b,200)    
resitev = []
for i in y:
    tren = []
    for j in x:
        tren.append(uiks(j,i,b))
    resitev.append(tren)
    
    
X, Y = np.meshgrid(x,y)
CS = plt.contourf(X,Y,resitev,50)
cbar = plt.colorbar(CS)
plt.title(r"Hitrostni profil, $\frac{a}{b} = 5$")
plt.hot()
plt.savefig("profilcek4.pdf")
plt.show()    
"""    
    
    
    
"""
x = np.linspace(1,2,300)
y = lenaloga(20,20,x)
yy = canal(x,300)
plt.plot(x,y,label="Razvoj po sinusih")
plt.plot(x,yy,label="Iz navodil")
plt.ylabel(r"C")
plt.xlabel(r"$\frac{a}{b}$")
plt.legend()
plt.savefig("moja5.pdf")
plt.show()
"""
"""
x = range(1,50)
plt.plot(x,[lenaloga(50,m,1) for m in x],label=r"$\frac{a}{b}=1$")
plt.plot(x,[lenaloga(50,m,2) for m in x],label=r"$\frac{a}{b}=2$")
plt.plot(x,[lenaloga(50,m,0.5) for m in x],label=r"$\frac{a}{b}=0.5$")
plt.plot(x,[lenaloga(50,m,5) for m in x],label=r"$\frac{a}{b}=5$")
plt.plot(x,[lenaloga(50,m,0.2) for m in x],label=r"$\frac{a}{b}=0.2$")
plt.title("C v odvisnosti od števila seštetih členov po n, m=200")
plt.legend()
plt.savefig("moja2.pdf")
plt.show()
"""



    
"""
vred = []
kot = 0
while kot<=3.15:
    temp = []
    radij=0
    while radij<=1.00002:
        temp.append(vrednost(radij,kot))
        radij+=0.01
    vred.append(temp)
    kot+=0.05
    
r = [0]
while r[-1]<1:
    r.append(r[-1]+0.01)
p = [0]
while p[-1]<np.pi:
   p.append(p[-1]+0.05)
R, P = np.meshgrid(r,p)
X,Y = R*np.cos(P), R*np.sin(P)

CS = plt.contourf(X,Y,vred,levels=50)
cbar = plt.colorbar(CS)
plt.title("Hitrostni profil")
plt.hot()
plt.savefig("profilkajwtf.pdf")
plt.show()
"""
    
    
    
    
    