# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import linalg as lin
from scipy.optimize import fsolve

#0 do 1 razdelimo na N+1 x0,x1,...xN točk
#velikost matrike je N-1 za točke x1,x2,...xN-1 torej mora bit N-1 točk




def analiticna(x):
    ksi1 = 0.37929
    ksi2= 2.73468
    temp = np.cosh(ksi2*(1-2*x))/np.cosh(ksi2)
    return -2*np.log(temp)


def newton(x0,N,st):
    A = itermatrika2(N)
    h = 1/N
    x = x0
    for i in range(st):
        fun = A*x - 0.5*h**2 *np.exp(x)
        fun = -1*fun
        jak = jakobi(x,N)
        delta = lin.solve_banded((1,1),jak,fun)
        x = x+delta
    return x

    
    
def itermatrika(w,N):
    h = 1/N
    matrika = []
    for i in range(N-1): #vrstica
        matrika.append([])
        for j in range(N-1): #stolpec
            if i==0: #prvo vrstico posebe
                if j==0:
                    matrika[i].append(w)
                elif j==1:
                    matrika[i].append(0.5)
                else:
                    matrika[i].append(0)
            elif i==N-2:
                if j==N-2:
                    matrika[i].append(w)
                elif j==N-3:
                    matrika[i].append(0.5)
                else:
                    matrika[i].append(0)                
            else:
                if i==j:
                    matrika[i].append(w)
                elif j==i-1 or j==i+1:
                    matrika[i].append(0.5)
                else:
                    matrika[i].append(0)
    return (1/(1+w)*np.matrix(matrika))

def itermatrika2(N):
    h = 1/N
    matrika = []
    for i in range(N-1): #vrstica
        matrika.append([])
        for j in range(N-1): #stolpec
            if i==0: #prvo vrstico posebe
                if j==0:
                    matrika[i].append(-2)
                elif j==1:
                    matrika[i].append(1)
                else:
                    matrika[i].append(0)
            elif i==N-2:
                if j==N-2:
                    matrika[i].append(-2)
                elif j==N-3:
                    matrika[i].append(1)
                else:
                    matrika[i].append(0)                
            else:
                if i==j:
                    matrika[i].append(-2)
                elif j==i-1 or j==i+1:
                    matrika[i].append(1)
                else:
                    matrika[i].append(0)
    return (-1*0.5)*np.matrix(matrika)    
"""    
def jakobi(x,N):
    h = 1/N
    matrika = []
    for i in range(N-1): #vrstica
        matrika.append([])
        for j in range(N-1): #stolpec
            if i==0: #prvo vrstico posebe
                if j==0:
                    matrika[i].append(1+0.5*h**2 * np.exp(float(x[0])))
                elif j==1:
                    matrika[i].append(-0.5 + 0.5*h**2 *np.exp(float(x[0])))
                else:
                    matrika[i].append(0)
            elif i==N-2:
                if j==N-2:
                    matrika[i].append(1+0.5*h**2 * np.exp(float(x[N-2])))
                elif j==N-3:
                    matrika[i].append(-0.5+0.5*h**2 * np.exp(float(x[N-2])))
                else:
                    matrika[i].append(0)                
            else:
                if i==j:
                    matrika[i].append(1+0.5*h**2 * np.exp(float(x[i])))
                elif j==i-1 or j==i+1:
                    matrika[i].append(-0.5+0.5*h**2 * np.exp(float(x[i])))
                else:
                    matrika[i].append(0)
    return np.matrix(matrika)       
"""    

def jakobi(x,N):
    h = 1/N
    matrika = [[0],[],[]]
    for i in range(N-2):
        matrika[0].append(-0.5+0.5*h**2*np.exp(float(x[i])))
    for i in range(N-2):
        matrika[2].append(-0.5+0.5*h**2 * np.exp(float(x[i])))
    matrika[2].append(0)
    for i in range(N-1):
        matrika[1].append(1+0.5*h**2 * np.exp(float(x[i])))
    return np.matrix(matrika)       
def iteracija(x,w,N,A,D): #vrne naslednji vektor
    h = 1/N
    return A*x + h**2 /(2+2*w)*D*np.exp(x)
    


def vrni(w,N,st,D):
    u0 = np.linspace(0,1,N+1)
    u0 = u0[1:-1] #prvi in zadnji element sta 0
    u0 = list(map(lambda x : 16*x*(1-x),u0))
    u0 = np.matrix(u0).T    
    A = itermatrika(w,N)
    for i in range(st):
        u0 = iteracija(u0,w,N,A,D)
    x = np.linspace(0,1,N+1)
    u0 = np.array(u0)
    u0 = np.append(u0,0)
    u0 = np.insert(u0,0,0)
    return [x,u0]

def vrni2(N,st):
    u0 = np.linspace(0,1,N+1)
    u0 = u0[1:-1] #prvi in zadnji element sta 0
    u0 = list(map(lambda x : 16*x*(1-x),u0))
    u0 = np.matrix(u0).T    
    u0 = newton(u0,N,st)
    x = np.linspace(0,1,N+1)
    u0 = np.array(u0)
    u0 = np.append(u0,0)
    u0 = np.insert(u0,0,0)
    return [x,u0]



x = np.linspace(0,1,200)

#y je oblike [y,p]
def funkc(x,y):
    return [y[1],-np.exp(y[0])]

def funkcc(x):
    r = ode(funkc).set_integrator("dopri5")
    r.set_initial_value([0,x])
    t=0.1
    while t<1:
        xx = r.integrate(t)
        t+=0.1
    return xx[0]
aa = fsolve(funkcc,10)

def resi(a):
    x = [0]
    y = [0]
    t=0
    r = ode(funkc).set_integrator("dopri5")
    r.set_initial_value([0,a])
    while x[-1]<1:
        t+=0.05
        x.append(t)
        y.append(r.integrate(t)[0])
    return [x,y]                



a = resi(aa)
plt.title("Absolutne napake za več metod")
plt.plot(a[0],abs(np.array(a[1])-list(map(analiticna,a[0]))),label="Strelska")
a = vrni(0,1000,10000,1)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="Iterativna shema")
a = vrni2(100,5)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="Newtnova")
plt.xlim(0,1)
plt.legend()
plt.yscale("log")
plt.savefig("strelskaprimer2.pdf")
plt.show()




"""
plt.title(r"Rešitev za več $\delta$")
a = vrni(0,1000,10000,0.5)
plt.plot(a[0],a[1],label=r"$\delta=0.5$")
a = vrni(0,1000,10000,1)
plt.plot(a[0],a[1],label=r"$\delta=1$")
a = vrni(0,1000,10000,1.5)
plt.plot(a[0],a[1],label=r"$\delta=1.5$")
a = vrni(0,1000,10000,2)
plt.plot(a[0],a[1],label=r"$\delta=2$")
a = vrni(0,1000,10000,2.5)
plt.plot(a[0],a[1],label=r"$\delta=2.5$")
plt.legend()
plt.savefig("delte2.pdf")
plt.show()
"""






#PLOTI PRI VEČ N NEWTon
"""
plt.title(r"Rešitev pri več velikosti matrik, Newton, Št. iteracij=5")
#a = vrni2(50,100)
#plt.plot(a[0],a[1],label="N=50")
a = vrni2(10,5)
plt.plot(a[0],a[1],label="N=10")
#a = vrni2(25,100)
#plt.plot(a[0],a[1],label="N=25")
a = vrni2(100,5)
plt.plot(a[0],a[1],label="N=100")
#a = vrni2(200,100)
#plt.plot(a[0],a[1],label="N=200")
plt.plot(x,list(map(analiticna,x)),label="Prava vrednost")
plt.legend(loc="upper right")
plt.savefig("18.pdf")
plt.show()
"""

"""
plt.title(r"Absolutne napake pri več velikosti matrik,Newton , 5 iteracij")
#a = vrni2(50,1000)
#plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="N=50")
a = vrni2(10,5)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="N=10")
#a = vrni2(25,1000)
#plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="N=25")
a = vrni2(100,5)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="N=100")
#a = vrni2(200,1000)
#plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="N=200")
#plt.plot(x,list(map(analiticna,x)),label="Prava vrednost")
plt.legend(loc="upper right")
plt.yscale("log")
plt.savefig("19.pdf")
plt.show()
"""








#PLOTI PRI VEČ ST NEWTON
"""
x = np.linspace(0,1,100)
plt.title(r"Rešitev pri več številih iteracij I,Newton, N=100")
#a = vrni2(100,500)
#plt.plot(a[0],a[1],label="I=500")
#a = vrni2(100,1000)
#plt.plot(a[0],a[1],label="I=1000")
a = vrni2(100,1)
plt.plot(a[0],a[1],label="I=1")
a = vrni2(100,5)
plt.plot(a[0],a[1],label="I=5")
a = vrni2(100,10)
plt.plot(a[0],a[1],label="I=10")
plt.plot(x,list(map(analiticna,x)),label="Prava vrednost")
plt.legend(loc="upper right")
plt.savefig("21.pdf")
plt.show()
"""


#NAPAKE PRI VEČ ST NEWTON
"""
plt.title(r"Absolutne napake pri več številih iteracij I, Newton, N=100")
#a = vrni2(100,500)
#plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=500")
#a = vrni2(100,1000)
#plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=1000")
#a = vrni2(100,5000)
#plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=5000")
a = vrni2(100,1)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=1")
a = vrni2(100,5)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=5")
a = vrni2(100,10)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=10")
#plt.plot(x,list(map(analiticna,x)),label="Prava vrednost")
plt.legend(loc="upper right")
plt.yscale("log")
plt.savefig("22.pdf")
plt.show()
"""






#ABS napakae pri več w
"""
plt.title(r"Absolutne napake pri $\omega$, I=500, N=1000")
a = vrni(0.5,1000,500)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label=r"$\omega = 0.5$")
a = vrni(0.8,1000,500)
#plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label=r"$\omega = 0.8$")
a = vrni(1,1000,500)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label=r"$\omega = 1$")
a = vrni(0.2,1000,500)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label=r"$\omega = 0.2$")
a = vrni(0,1000,500)
#plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label=r"$\omega = 0$")

#plt.plot(x,list(map(analiticna,x)),label="Prava vrednost")
plt.legend(loc="upper right")
plt.yscale("log")
#plt.savefig("17.pdf")
plt.show()
"""

#PLOTI PRI VEČ ŠT
"""
plt.title(r"Rešitev robnega problema pri več številih iteracij I, $\omega = 0.5$, N=1000")
a = vrni(0.5,1000,500)
plt.plot(a[0],a[1],label="I=500")
a = vrni(0.5,1000,1000)
plt.plot(a[0],a[1],label="I=1000")
a = vrni(0.5,1000,5000)
plt.plot(a[0],a[1],label="I=5000")
a = vrni(0.5,1000,10000)
plt.plot(a[0],a[1],label="I=10000")
a = vrni(0.5,1000,20000)
plt.plot(a[0],a[1],label="I=20000")
plt.plot(x,list(map(analiticna,x)),label="Prava vrednost")
plt.legend(loc="upper right")
plt.savefig("15.pdf")
plt.show()
"""



#ABS NAP PRI VEČ ŠT
"""
plt.title(r"Absolutne napake pri več številih iteracij I, $\omega = 0.5$, N=50")
a = vrni(0.5,1000,500)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=500")
a = vrni(0.5,1000,1000)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=1000")
a = vrni(0.5,1000,5000)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=5000")
a = vrni(0.5,1000,10000)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=10000")
a = vrni(0.5,1000,20000)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="I=20000")
#plt.plot(x,list(map(analiticna,x)),label="Prava vrednost")
plt.legend(loc="upper right")
plt.yscale("log")
plt.savefig("16.pdf")
plt.show()
"""











#PLOTI PRI VEČ N
"""
plt.title(r"Rešitev robnega problema, $\omega = 0.5$, 10000 iteracij")
a = vrni(0.5,500,10000)
plt.plot(a[0],a[1],label="N=500")
a = vrni(0.5,100,10000)
plt.plot(a[0],a[1],label="N=100")
a = vrni(0.5,50,10000)
plt.plot(a[0],a[1],label="N=50")
a = vrni(0.5,10,10000)
plt.plot(a[0],a[1],label="N=10")
a = vrni(0.5,1000,10000)
plt.plot(a[0],a[1],label="N=1000")
plt.plot(x,list(map(analiticna,x)),label="Prava vrednost")
plt.legend(loc="upper right")
plt.savefig("13.pdf")
plt.show()
"""

#ABSOLUTNE NAPAKE PRI VEČ N
"""
plt.title(r"Absolutne napake, $\omega = 0.5$, 10000 iteracij")
a = vrni(0.5,500,10000)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="N=500")
a = vrni(0.5,100,10000)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="N=100")
a = vrni(0.5,50,10000)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="N=50")
a = vrni(0.5,10,10000)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="N=10")
a = vrni(0.5,1000,10000)
plt.plot(a[0],abs(a[1]-list(map(analiticna,a[0]))),label="N=1000")
#plt.plot(x,list(map(analiticna,x)),label="Prava vrednost")
plt.legend(loc="upper right")
plt.yscale("log")
plt.savefig("14.pdf")
plt.show()
"""















#PREVERJA KOLIKO HITRO KONVERGIRA
"""
u0 = np.linspace(0,1,51)
u0 = u0[1:-1] #prvi in zadnji element sta 0
xxx = list(map(analiticna,u0))
u0 = list(map(lambda x : x*(1-x),u0))
u0 = np.matrix(u0).T    
A = itermatrika(1,50)
for i in range(30000):
    u0 = iteracija(u0,1,50,A)
    if abs(sum(u0)-sum(xxx))<0.01:
        print(i)
        break
"""
#NARIŠE TO
"""
plt.title(r"Število iteracij za abs. natančnost 0.01, N=50")
plt.ylabel("Št iteracij")
plt.xlabel(r"$\omega$")
plt.plot([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],[3383,3721,4060,4398,4737,5076,5414,5753,6091,6430,6769],"o-")
plt.savefig("8.pdf")
plt.show()
"""
