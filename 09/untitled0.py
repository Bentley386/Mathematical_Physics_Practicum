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

D = 0.000001

def zac(x):
    return 20*np.exp(-(x-0.05)**2/0.01)
def diskretiziraj(L,N):
    x = []
    rejt = 1/(N)*L
    for i in range(N):
        x.append(rejt*i)
    return x
def four(T,x):
    N = len(x)
    koef = []
    for i in range(N):
        suma=0
        for j in range(N):
            suma += T[j]*np.exp(2*np.pi*1j*x[j]*i) 
        koef.append(suma)
    return [x,koef]
def ifour(T,x):
    N = len(x)
    koef = []
    for i in range(N):
        suma=0
        for j in range(N):
            suma += T[j]*np.exp(-2*np.pi*1j*x[j]*i) 
        koef.append(1/N*suma)
    return [x,koef]    
def euler(T0,k,h):
    T = [T0]
    t = [0]
    cas = 0
    i = 0
    temp = T0
    while cas < 100:
        cas+=h
        i+=1
        temp = temp -h*D*4*(np.pi**2 * k**2)*temp
        if i==10:
            i=0
            T.append(temp)
            t.append(cas)
    return [cas,temp]

def euler2pomozna(T0,k,h,N,t):
    cas = 0
    temp = T0
    while cas < t:
        cas+=h
        temp = temp -h*D*4*(np.pi**2 * (k/N)**2)*temp
    return temp
    
def euler2(T0,t,h):
    koef = []
    N = len(T0)
    for i in range(len(T0)):
        koef.append(euler2pomozna(T0[i],i,h,N,t))
    return koef

def anal(zacetna,t):
    N = len(zacetna)
    resitve =[]
    for i in range(len(zacetna)):
        resitve.append(zacetna[i]*np.exp(-4*D*np.pi**2*(i/N)**2*t))
    return resitve
    
def naredi(N,t,h):
    x = diskretiziraj(0.1,N)
    zacetna = list(map(zac,x))
    zacetnaTrans = np.fft.fft(zacetna)
    #temperature = euler2(zacetnaTrans,t,h)
    temperature=anal(zacetnaTrans,t)
    x.append(0.1)
    rezul=np.real(np.fft.ifft(temperature))
    rezul = np.append(rezul,rezul[0])
    return [x,rezul]

def naredi2(N,t,h):
    x = diskretiziraj(0.1,N)
    x.append(0.1)
    tt = [i for i in range(t)]
    vrednosti = []
    for i in tt:
        vrednosti.append(naredi(N,i,h)[1])
    return [x,tt,vrednosti]

#v prvi vrstici imamo od leve proti desni pri x od 0 do 1,
#po stolpcih dol pol cas od 0 do 100

       
"""
a = naredi2(20,1000,1)
x,t = np.meshgrid(a[0],a[1])

fig = plt.figure()
ax = fig.gca(projection="3d")
plt.xlim(0,0.1)
plt.xlabel("x")
plt.ylabel("t")
ax.plot_surface(x,t,a[2],rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0)
plt.savefig("3d4.pdf")
plt.show()
"""

"""
plt.title("T(x) pri več časih, h=1,N=20")
plt.xlim(0,0.1)
a = naredi(20,0,1)    
plt.plot(a[0],a[1],label="t=0")
a = naredi(20,20,1)    
plt.plot(a[0],a[1],label="t=20")
a = naredi(20,100,1)    
plt.plot(a[0],a[1],label="t=100")
a = naredi(20,500,1)    
plt.plot(a[0],a[1],label="t=500")
a = naredi(20,1000,1)    
plt.plot(a[0],a[1],label="t=1000")
a = naredi(20,2000,1)    
plt.plot(a[0],a[1],label="t=2000")
plt.legend()
plt.savefig("veccasih.pdf")
plt.show()
"""








def zacetni(x):
    return np.sin(np.pi*x)
x = diskretiziraj(1,20)
x.append(1) #dolzina zdaj 21
zac = list(map(zacetni,x))
zac=zac[1:-1]
def bspline(xk,x,h):
    if x<= (xk-2*h) or x>= (xk+2*h):
        return 0
    elif x<= (xk-h) and x>= (xk-2*h):
        return 1/(6*h**3)*(x-(xk-2*h))**3
    elif x<= (xk+2*h) and x>= (xk+h):
        return -1/(6*h**3)*(x-(xk+2*h))**3
    elif x<= xk and x>= (x-h):
        return 1/6+1/(2*h)*(x-(xk-h))+1/(2*h**2)*(x-(x-h))**2 -1/(2*h**3)*(x-(x-h))**3
    elif x>= xk and x<= (x+h):
        return 1/6-1/(2*h)*(x-(xk+h))+1/(2*h**2)*(x-(x+h))**2 +1/(2*h**3)*(x-(x+h))**3

A=[]
for i in range(1,20):
    vrstica = []
    for j in range(1,20):
        if i==1:
            if j==1:
                vrstica.append(4)
            elif j==2:
                vrstica.append(1)
            else:
                vrstica.append(0)
        elif i==19:
            if j==19:
                vrstica.append(4)
            elif j==18:
                vrstica.append(1)
            else:
                vrstica.append(0)
        else:
            if i==j:
                vrstica.append(4)
            elif abs(i-j)==1:
                vrstica.append(1)
            else:
                vrstica.append(0)
    A.append(vrstica)
A = np.matrix(A)
        
B=[]
for i in range(1,20):
    vrstica = []
    for j in range(1,20):
        if i==1:
            if j==1:
                vrstica.append(-2)
            elif j==2:
                vrstica.append(1)
            else:
                vrstica.append(0)
        elif i==19:
            if j==19:
                vrstica.append(-2)
            elif j==18:
                vrstica.append(1)
            else:
                vrstica.append(0)
        else:
            if i==j:
                vrstica.append(-2)
            elif abs(i-j)==1:
                vrstica.append(1)
            else:
                vrstica.append(0)
    B.append(vrstica)
B = 6*D/((1/20)**2)*np.matrix(B)
 
a0 = solve(A,6*(np.matrix(zac).T))

cas = 0
t=[0]
aji=[a0]
while cas<5000:
    cas+=1
    t.append(cas)
    novi = solve(A-0.5*B,(A+0.5*B)*aji[-1])
    aji.append(novi)
    
for i in range(len(aji)):
    aji[i] = np.append(aji[i],0)
    aji[i] = np.append(aji[i],aji[i][-2])
    aji[i] = np.insert(aji[i],0,0)
    aji[i] = np.insert(aji[i],0,aji[i][1])    

def temper(iks,t):
    suma = 0
    os = x[:]
    os.insert(0,-1/20)
    os.append(1+1/20)
    for i in range(len(os)):
        suma+=t[i]*bspline(os[i],iks,1/20)
    return suma
    #x je diskretiziraj 1,20 append 1
T = []
for i in aji:
    vrstica = []
    for j in x:
        vrstica.append(temper(j,i))
    T.append(vrstica)
    
x,t = np.meshgrid(x,t)

fig = plt.figure()
ax = fig.gca(projection="3d")
plt.title("Metoda z B-zlepki, D=0.000001,L=1")
plt.xlabel("x")
plt.ylabel("t")
ax.plot_surface(x,t,T,rstride=1,cstride=1,cmap=cm.coolwarm,linewidth=0)
plt.savefig("splajn4.pdf")
#plt.show()