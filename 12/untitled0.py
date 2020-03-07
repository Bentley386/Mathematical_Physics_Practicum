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
from scipy.special import beta
import scipy.sparse
from scipy.integrate import quad
from scipy.integrate import complex_ode

def A(mm,nn,m,n):
    if m!=mm:
        return 0
    return -np.pi/2 *n*nn*(3+4*m)/(2+4*m+n+nn)*beta(n+nn-1,3+4*m)
def b(m,n):
    return -2/(2*m+1) * beta(2*m+3,n+1)
    
    
def sestavi(mm,nn,m,n):
    matrike = []
    for i in range(mm):
        matrika = []
        for j in range(1,nn+1):
            vrstica = []
            for k in range(1,n+1):
                vrstica.append(A(i,j,i,k))
            matrika.append(vrstica)
        matrike.append(np.matrix(matrika))
    return scipy.sparse.block_diag(matrike,"csr")
            
def sestavib(mm,nn):
    vektor = []
    for m in range(mm):
        for n in range(1,nn+1):
            vektor.append(b(m,n))
    return np.matrix(vektor).T


def c(m,n):
    matrika = sestavi(m,n,m,n)
    vektor = sestavib(m,n)
    x = scipy.sparse.linalg.spsolve(matrika,vektor)
    x = np.matrix(x).T
    return -32/np.pi * (vektor.T * x)
print(c(50,50))
#matrika = sestavi(10,10,10,10)
#vektor = sestavib(10,10)
#a = scipy.sparse.linalg.spsolve(matrika.T,vektor)
#print(a)
"""
def vrednost(a,x,fi):
    i=0
    suma = 0
    for m in range(10):
        for n in range(1,11):
            clen= a[i]*np.sin((2*m+1)*fi)*x**(2*m+1)*(1-x)**n
            suma+=clen
            i+=1
    return suma

            
    
plt.title("C pri več seštetih členih")    
indeksi = range(1,20)
cji = [float(c(1,i)) for i in indeksi]
plt.plot(indeksi,cji,label="po m=1")
indeksi = range(1,20)
cji = [float(c(2,i)) for i in indeksi]
plt.plot(indeksi,cji,label="po m=2")
indeksi = range(1,20)
cji = [float(c(3,i)) for i in indeksi]
plt.xlabel("Število seštetih členov po n")
plt.ylabel("C")
plt.plot(indeksi,cji,label="po m=3")
indeksi = range(1,20)
cji = [float(c(4,i)) for i in indeksi]
plt.plot(indeksi,cji,label="po m=4")
indeksi = range(1,20)
cji = [float(c(5,i)) for i in indeksi]
plt.plot(indeksi,cji,label="po m=5")
plt.legend(loc="lower right")
plt.savefig("2.pdf")
plt.show()
"""

"""
vred = []
kot = 0
while kot<=3.15:
    temp = []
    radij=0
    while radij<=1.00002:
        temp.append(vrednost(a,radij,kot))
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

CS = plt.contourf(X,Y,vred,50)
cbar = plt.colorbar(CS)
plt.title("Hitrostni profil")
plt.hot()
plt.savefig("profilkaj.pdf")
plt.show()


def analkoef(k,t):
    return np.sin(k*np.pi/2) * jv(k,np.pi) * np.exp(1j*k*t)
def anal(x,t):
    return np.sin(np.pi * np.cos(x+t))
  
def zacetni(k):
    def pomozna(x,k):
        return 1/(2*np.pi)* np.sin(np.pi*np.cos(x))*np.cos(k*x)
    def pomozna2(x,k):
        return 1/(2*np.pi)*np.sin(np.pi*np.cos(x))*np.sin(k*x)
    return quad(pomozna,0,2*np.pi,k)[0] +1j * quad(pomozna2,0,2*np.pi,k)[0]
    #return quad(pomozna,0,2*np.pi,k)[0]

def razvoj(k,t):
    def odv(t,y):
        return 1j*k*y
    integ = complex_ode(odv)
    integ.set_integrator("dopri5")
    integ.set_initial_value(zacetni(k))
    t = list(range(0,t+1))
    y = []
    for i in t:
        y.append(integ.integrate(i))
    return [t,y]

def aj(k,t):
    def odv(t,y):
        return 1j*k*y
    integ = complex_ode(odv)
    integ.set_integrator("dopri5",nsteps=1000)
    integ.set_initial_value(zacetni(k))
    return integ.integrate(t)

def u(x,t,N):
    i = -N/2
    suma = 0
    while i<=N/2:
        suma+= aj(i,t)*np.exp(1j*i*x)
        i+=1
    return suma

    
    
    
    
    

plt.title("abs. napake za več N za t=15")
x = np.linspace(0,2*np.pi)
y = abs(anal(x,15)-u(x,15,6))
plt.plot(x,y,label="N=6")
y = abs(anal(x,15)-u(x,15,10))
plt.plot(x,y,label="N=10")
y = abs(anal(x,15)-u(x,15,14))
plt.plot(x,y,label="N=14")
y = abs(anal(x,15)-u(x,15,20))
plt.plot(x,y,label="N=20")
y = abs(anal(x,15)-u(x,15,100))
plt.plot(x,y,label="N=100")
plt.xlabel("x")
plt.yscale("log")
plt.legend()
#plt.savefig("napaka3.pdf")
plt.show()
"""





    
"""
aproks = razvoj(-1,50)
realno = [list(range(0,51)),[analkoef(-1,i) for i in range(0,51)]]
plt.title(r"Koeficient $a_{-1} (t)$")
plt.xlabel("t")
plt.plot(aproks[0],aproks[1],label="aproksimacija")
plt.plot(realno[0],realno[1],label="analitično")
plt.legend()
plt.savefig("koef/minus1.pdf")
plt.show()


def sestavi(n,k):
    h = 2*np.pi/(n-1)
    matrika = []
    for i in range(n):
        vrstica = []
        for j in range(n):
            if i==j:
                vrstica.append(1-k/h)
            elif j == i+1:
                vrstica.append(k/h)
            elif i==n-1 and j==0:
                vrstica.append(k/h)
            else:
                vrstica.append(0)
        matrika.append(vrstica)
    return np.matrix(matrika)
    
def initial(x):
    return np.sin(np.pi*np.cos(x))    
    

def casrazvoj(t,y0,n,k):
    matrika = sestavi(n,k)
    t0 = 0
    y0 = np.matrix(y0).T
    while t0<=t:
        t0+=k
        y0 = matrika*y0
    return y0.T

plt.title("Absolutna razlika difrenčne metode in metode Galerkina, N=20,n=10")    
x = np.linspace(0,2*np.pi,10)
y0 = initial(x)
y0 = casrazvoj(5,y0,10,0.01)
prava = u(x,5,20)
plt.plot(x,abs(np.array(y0[0]).flatten()-prava),label="m=10")

x = np.linspace(0,2*np.pi,10)
y0 = initial(x)
y0 = casrazvoj(5,y0,10,0.1)
prava = u(x,5,20)
plt.plot(x,abs(np.array(y0[0]).flatten()-prava),label="m=50")
x = np.linspace(0,2*np.pi,10)
y0 = initial(x)
y0 = casrazvoj(5,y0,10,0.2)
prava = u(x,5,20)
plt.plot(x,abs(np.array(y0[0]).flatten()-prava),label="m=25")

plt.legend()
plt.yscale("log")
#plt.savefig("nova3.pdf")
plt.show()
    
    
"""
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    