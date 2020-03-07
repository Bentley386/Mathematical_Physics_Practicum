# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt


#x vektor [x,p]
def funk1(t,x):
    return [x[1],-np.sin(x[0])]

def dod2(t,x):
    return[x[1],10*np.cos(t)-x[0]+x[1]*(1-x[0]**2)]
def dod(t,x,w):
    return [x[1],1.5*np.cos(w*t)-np.sin(x[0]) - 0.5*x[1]]
t = np.linspace(0,10,100)
def runge(x0,u0,t,h):
    tt = [0]
    while tt[-1]<t:
        tt.append(tt[-1]+h)
    r = ode(dod2).set_integrator("dopri5")
    r.set_initial_value([x0,u0])
    #r.set_f_params(0.05)
    y = []
    for i in tt:
        y.append(r.integrate(i)[0])
    return [tt,y]

def runge3(w):
    tt = [0]
    while tt[-1]<100:
        tt.append(tt[-1]+0.1)
    r = ode(dod).set_integrator("dopri5")
    r.set_f_params(w)
    r.set_initial_value([1,0])
    y = []
    for i in tt:
        y.append(r.integrate(i)[0])
    return [tt,y]
def runge2(x0,u0,t,h):
    tt = [0]
    while tt[-1]<t:
        tt.append(tt[-1]+h)
    r = ode(dod2).set_integrator("dopri5")
    r.set_initial_value([x0,u0])
    y = []
    for i in tt:
        y.append(r.integrate(i)[1])
    return [tt,y]
    
def trapezna(x0,u0,t,h):
    x5 = ode(funk1).set_integrator("dopri5").set_initial_value([x0,u0]).integrate(h)[0]
    u5 = ode(funk1).set_integrator("dopri5").set_initial_value([x0,u0]).integrate(h/2)[1]
    tx = [0,h]
    x = [x0,x5]
    tu = [0,h/2]
    u = [u0,u5]
    while tx[-1] < t:
        tu.append(tu[-1]+h)
        u.append(u[-1]-h*np.sin(x[-1]))
        tx.append(tx[-1]+h)
        x.append(x[-1] + h* u[-1])
    return [tx,x,tu,u]

def trapezna2(x0,u0,t,h):
    x5 = ode(funk1).set_integrator("dopri5").set_initial_value([x0,u0]).integrate(3*h/2)[0]
    u5 = ode(funk1).set_integrator("dopri5").set_initial_value([x0,u0]).integrate(h)[1]
    tx = [0,3*h/2]
    x = [x0,x5]
    tu = [0,h]
    u = [u0,u5]
    while tu[-1] < t:
        tu.append(tu[-1]+h)
        u.append(u[-1]-h*np.sin(x[-1]))
        tx.append(tx[-1]+h)
        x.append(x[-1] + h* u[-1])
    return [tx,x,tu,u]

def abs_napaka(x0,u0,t,h):
    prvaS = runge(x0,u0,t,h)
    drugaS = trapeznaa(x0,u0,t,h)
    return [prvaS[0],abs(np.array(prvaS[1])-np.array(drugaS[1]))]
def abs_napaka2(x0,u0,t,h):
    prvaS = runge2(x0,u0,t,h)
    drugaS = trapeznaa2(x0,u0,t,h)
    return [prvaS[0],abs(np.array(prvaS[1])-np.array(drugaS[3]))]


def trapeznaa(x0,u0,t,h):
    x5 = ode(funk1).set_integrator("dopri5").set_initial_value([x0,u0]).integrate(h)[0]
    u5 = ode(funk1).set_integrator("dopri5").set_initial_value([x0,u0]).integrate(h/2)[1]
    tx = [0,h]
    x = [x0,x5]
    tu = [0,h/2]
    u = [u0,u5]
    while tx[-1] < t:
        tu.append(tu[-1]+h)
        tx.append(tx[-1]+h)
        tempu=0.5
        for i in range(4):
            tempu=u[-1]-h*np.sin(x[-1])+h/12*(2*np.sin(x[-1])-np.sin(x[-2])-np.sin(x[-1]+h*tempu))
        u.append(tempu)
        x.append(x[-1]+h*u[-1])
    return [tx,x,tu,u]    

def trapeznaa2(x0,u0,t,h):
    x55 = ode(funk1).set_integrator("dopri5").set_initial_value([x0,u0]).integrate(h/2)[0]
    x5 = ode(funk1).set_integrator("dopri5").set_initial_value([x0,u0]).integrate(3*h/2)[0]
    u5 = ode(funk1).set_integrator("dopri5").set_initial_value([x0,u0]).integrate(h)[1]
    tx = [h/2,3*h/2]
    x = [x55,x5]
    tu = [0,h]
    u = [u0,u5]
    while tu[-1] < t:
        tu.append(tu[-1]+h)
        tx.append(tx[-1]+h)
        tempu=0.5
        for i in range(4):
            tempu=u[-1]-h*np.sin(x[-1])+h/12*(2*np.sin(x[-1])-np.sin(x[-2])-np.sin(x[-1]+h*tempu))
        u.append(tempu)
        x.append(x[-1]+h*u[-1])
    return [tx,x,tu,u]    

frekvenca = [0.01]
ips = [max(runge3(0.01)[1])]
while frekvenca[-1] < 2:
    ips.append(max(runge3(frekvenca[-1])[1]))
    frekvenca.append(frekvenca[-1]+0.01)

"""
plt.plot(frekvenca,ips)
plt.title("Resonančna krivulja , v=1.5, x(0)=1,v(0)=0")
plt.xlabel(r"Frekvenca $\omega_0$")
plt.ylabel("Amplituda  oscilacije")
plt.savefig("reso/6.pdf")    
plt.show()
"""

#plt.plot(t,y,label="Prava vrednost")


"""
x = abs_napaka2(1,1,150,0.1)
plt.plot(x[0],x[1],label="h=0.1")
x = abs_napaka2(1,1,150,0.05)
plt.plot(x[0],x[1],label="h=0.05")
x = abs_napaka2(1,1,150,0.02)
plt.plot(x[0],x[1],label="h=0.02")
x = abs_napaka2(1,1,150,0.03)
plt.plot(x[0],x[1],label="h=0.03")
plt.yscale("log")
plt.legend(loc="lower right")
plt.title("Matematično nihalo - Trapezna shema v2, Abs. napake: x(0)=1,v(0)=1")
plt.ylabel("Napaka Hitrosti")
plt.xlabel("Čas")
plt.savefig("mat/33.pdf")
plt.show()
"""



x = runge(1,1,60,0.01)[1]
y = runge2(1,1,60,0.01)[1]
plt.xlabel("Odmik")
plt.ylabel("Hitrost")
plt.title("Fazni portret van der Pol: x(0)=1,v(0)=1,h=0.01, 10 nihajev")
plt.plot(x,y)
plt.savefig("pol6.pdf")
plt.show()







"""
prava = runge(1,1,50,0.01)
plt.plot(prava[0],prava[1],label="Prava vrednost")
#x = trapeznaa(1,3,40,1)
#plt.plot(x[0],x[1],label="h=1")
#x = trapeznaa(1,1,40,0.5)
#plt.plot(x[0],x[1],label="h=0.5")
#x = trapeznaa(1,1,40,0.1)
#plt.plot(x[0],x[1],label="h=0.1")
#plt.legend()
plt.title("Van der Polov oscilator- RK h=0.01: x(0)=1, v(0)=1")
plt.ylabel("Odmik")
plt.xlabel("Čas")
plt.savefig("pol4.pdf")
plt.show()
"""