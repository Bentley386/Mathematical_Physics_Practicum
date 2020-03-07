# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
from scipy.special import airy
from scipy.optimize import newton
import numpy as np
from scipy.optimize import fsolve

alfa = 0.355028053887817239
beta = 0.258819403792806798

def ai(x):
    return airy(x)[0]
def aio(x):
    return airy(x)[1]
def bi(x):
    return airy(x)[2]
def bio(x):
    return airy(x)[3]

def f(x,n=1000000,e=10**-12):
    a0=1
    y = 1
    for k in range(1,n):
        clen=a0*3*(x**3)*(k-2/3)/(3*k*(3*k-1)*(3*k-2))
        y+= clen
        if e!=0 and abs(clen)<=e :
            return y
            break   
        a0 = clen
    if e==0:
        return y
    if e!=0:
        print("vec kot miljon ponovitev")
    
def g(x,n=1000000,e=10**-12):
    a0=x
    y = x
    for k in range(1,n):
        clen=a0*3*(x**3)*(k-1/3)/(3*k*(3*k-1)*(3*k+1))
        y+= clen
        if e!= 0 and abs(clen)<=e:
            return y
            break   
        a0 = clen
    if e==0:
        return y
    if e!=0:
        print("vec kot miljon ponovitev")
        
def L(x,n=1000000,e=10**-12):
    a0=1
    y = 1
    for k in range(1,n):
        clen=a0*(3*k-2.5)*(3*k-1.5)*(3*k-0.5)/(54*x*k*(k-0.5))
        y+= clen
        if e!= 0 and abs(clen)<=e:
            return y
            break
        if abs(clen)>abs(a0):
            return y
        a0 = clen
    if e==0:
        return y
    if e!=0:
        print("vec kot miljon ponovitev")
        return y

def P(x,n=1000000,e=10**-10):
    a0=1
    y =1
    for k in range(1,n):
        clen=-a0*(6*k-5.5)*(6*k-4.5)*(6*k-3.5)*(6*k-2.5)*(6*k-1.5)*(6*k-0.5)\
        /(54**2*x**2*2*k*(2*k-1)*(2*k-1.5)*(2*k-0.5))
        y+= clen
        if e!= 0 and abs(clen)<=e:
            return y
            break
        if abs(clen)>abs(a0):
            return y
        a0 = clen
    if e==0:
        return y
    if e!=0:
        print("vec kot miljon ponovitev")
        return y
        
def Q(x,n=1000000,e=10**-10):
    a0=5/(72*x)
    y =5/(72*x)
    for k in range(1,n):
        clen=-a0*(6*k-2.5)*(6*k-1.5)*(6*k-0.5)*(6*k+0.5)*(6*k+1.5)*(6*k+2.5)\
        /(54**2*x**2*2*k*(2*k+1)*(2*k+0.5)*(2*k-0.5))
        y+= clen
        if e!= 0 and abs(clen)<=e:
            return y
            break
        if abs(clen)>abs(a0):
            return y
        a0 = clen
    if e==0:
        return y
    if e!=0:
        print("vec kot miljon ponovitev")
        return y
        
#Ai za majhne argumente    
def aim(x,n=1000000,e=10**-12):
    if x==0:
        return alfa
    return alfa*f(x,n,e)-beta*g(x,n,e)

#Bi za majhne argumente    
def bim(x,n=1000000,e=10**-12):
    if x==0:
        return alfa*(3**0.5)
    return (3**0.5)*(alfa*f(x,n,e)+beta*g(x,n,e))

#bi za velike argumetne
def biv(x,n=1000000,e=10**-12):
    ceta = 2/3*(abs(x)**(3/2))
    if x>0:
        A = np.exp(ceta)/(np.pi**0.5*(x**0.25))
        return A*L(ceta,n,e)
    elif x<0:
        A = 1/(np.pi**0.5*((-x)**0.25))
        return A*(-np.sin(ceta-np.pi/4)*P(ceta,n,e) + np.cos(ceta-np.pi/4)*Q(ceta,n,e))

#ai za velike argumente
def aiv(x,n=1000000,e=10**-12):
    ceta = 2/3*(abs(x)**(3/2))
    if x>0:
        A = np.exp(-ceta)/(2*(np.pi**0.5)*(x**0.25))
        return A*L(-ceta,n,e)
    elif x<0:
        A = 1/(np.pi**0.5*((-x)**0.25))
        return A*(np.sin(ceta-np.pi/4)*Q(ceta,n,e) + np.cos(ceta-np.pi/4)*P(ceta,n,e))

#airy
def mojai(x,n=0,e=0):
    if abs(x)<6:
        return aim(x,10000,10**-14)
    else:
        return aiv(x,3000,0)

#bairy
def mojbi(x,n=0,e=0):
    if abs(x)<7:
        return bim(x,3000,10**-12)
    else:
        return biv(x,100000,0)
       
    
#vrne seznam 
def rezultat(funk,x,n=10000,e=10**-12):
    y = []
    for i in list(x):
        y.append(funk(i,n,e))
    return y

    
def f(z):
    return z**(2/3)*(1+5/48*(z**(-2))-5/36*(z**(-4))+77125/82944*(z**(-6))-108056875/6967296*(z**(-8)))    
    
    
def prvih100a2():
    nicle = []
    for s in range(1,101):
        arg = 3*np.pi*(4*s-1)/8
        nicle.append(-f(arg))
    return nicle
def prvih100b2():
    nicle = []
    for s in range(1,101):
        arg = 3*np.pi*(4*s-3)/8
        nicle.append(-f(arg))
    return nicle

    
    


#x = np.linspace(-90,-95,100)

#plt.plot(x,ai(x))
#plt.xscale("log")
#plt.yscale("log")
#plt.show()


def dx(f, x):
    return abs(0-f(x))
 
def newtons_method(f, df, x0, e):
    delta = dx(f, x0)
    while delta > e:
        x0 = x0 - f(x0)/df(x0)
        delta = dx(f, x0)
    #print('Root is at:{}'.format(x0))
    #print('f(x) at root is:{}'.format(f(x0)))
    return x0

def pomozna(nicla,nicle):
    for i in nicle:
        if abs(nicla-i)<10**-1:
            return False
    return True
def prvih100a(step=0.1):
    i = -2.43
    nicle = [-2.336141997554901]
    while len(nicle)<100:
        nicla = newtons_method(ai,aio,i,10**-2)
        print(nicla)
        if nicla not in nicle and pomozna(nicla,nicle) and nicla<0:
            nicle.append(nicla)
        i-=step
    return nicle
    
def prvih100b(step=0.1):
    i = -1.17371
    nicle = [-1.17371]
    while len(nicle)<101:
        nicla = newtons_method(bi,bio,i,10**-2)
        print(nicla)
        if nicla not in nicle and pomozna(nicla,nicle) and nicla<0:
            nicle.append(nicla)
        i-=step
    return nicle
#DODATNA NALOGA=====================================
mojenicle = sorted(prvih100b(),reverse=True)[:100]
pravenicle = prvih100b2()

for i in range(100):
    print("formula:{}, Newton{}\n".format(pravenicle[i],mojenicle[i]))

"""
mojenicle = sorted(prvih100a(),reverse=True)
pravenicle = prvih100a2()

for i in range(100):
    print("formula:{}, Newton{}\n".format(pravenicle[i],mojenicle[i]))
"""


