# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import math as m

def exata(x):
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = (x[i] ** 2) / 2
    return y


"""
/*
    Método Lax–Friedrichs - LxF
*/
"""
def laxFriedrichs(uu,Nx,Nt,h,dt):
    global K
    u = np.zeros([Nx+1,Nt+1])
    u[:,:] = uu[:,:]
    for t in range(1,Nt):
        for i in range(1,Nx-1):
            u[i,t+1] = 0.5*(u[i+1,t]+u[i-1,t]) - (K*dt/(2*h))*(u[i+1,t]-u[i-1,t])
        u[0,t+1] = u[1,t+1]
        u[Nx-1,t+1] = u[Nx-2,t+1]
    return u

"""
/*
    Método Kurganov-Tadmor - KT
*/
"""
def KT(uu, Nx, Nt, dx, dt):
    u = np.zeros([Nx+1,Nt+1])
    u[:,:] = uu[:,:]
    for t in range(1, Nt):
        for j in range(1, Nx):
            u[j,t] = u[j,0] + dt * (-(H(u,j,t,dx) - H(u, j-1, t,dx))/dx)

    return u

def f(u):
    return m.sin(u)

def dfdu(u):
    return m.cos(u)

def a(u, j, t, dx):
    a1 = dfdu(uPlus(u, j, t, dx))
    a2 = dfdu(uMinus(u, j, t, dx))
    if a1 > a2:
        return a1
    return a1

def H(u, j, t, dx):
    return 0.5*(f(uPlus(u, j, t, dx)+uMinus(u, j, t, dx)) - a(u, j, t, dx)*(uPlus(u, j, t, dx) - uMinus(u, j, t, dx)))

def uPlus(u, j, t, dx):
    return u[j+1,t] - 0.5 * dx * u[j+1,t]

def uMinus(u, j, t, dx):
    return u[j,t] + 0.5 * dx * u[j,t]    

############################################################

def condicaoInicial(Nx,Nt,h):
    u = np.zeros([Nx+1,Nt+1])
    for i in range(0,Nx+1):
        u[i,0] = m.sin(i*h)
    return u


def condicoesContorno(u,Nx,Nt):
    #for t in range(0, Nt+1):
    #u[0,0] = m.sin(0*h)
    #u[Nx-1,0] = m.sin((Nx-1)*h)
    return u

#Entrada
A = 0 #Limite inferior para x
B = m.pi*2 #Limite superior para x
Nx = 102 #Quantidade de elementos + (inicio + final) -> 2
h = (B - A)/np.float(Nx) #Discretização no espaço
T = 4 #Tempo final
dt = 0.5*h #Discretização no tempo
Nt = np.int(T/dt) #Quantidade de iterações no tempo
x = np.linspace(A,B,Nx+1) #Para plot das aproximações
K = 1 #Coeficiente convectivo
xe = np.linspace(A,B,1000) #Para plot da solução exata
ye = exata(xe) #Para plot da solução exata

u = condicaoInicial(Nx,Nt,h)

u = condicoesContorno(u,Nx,Nt)

#u_laxFri = laxFriedrichs(u,Nx,Nt,h,dt)
u_kt = KT(u,Nx,Nt,h,dt)

tj = 0 #Instante de tempo desejado
t = int(tj/dt) #Índice correspondente ///
t = 1

#Plota o grafico
#plt.title('Exata x Aproximada (400 elementos, dt = 1.5h)')
plt.title('Exata x Aproximada ('+ str(Nx) +  ' elementos , dt = ' + str(dt/h) + 'h)')
#plt.xlim (-1 ,1)
#plt.ylim (0 ,1.1)
plt.grid()

#plt.plot(x,u_laxFri[:,t],'b-',label = 'Lax-Friedrichs')
plt.plot(x,u_kt[:,t],'b-',label = 'Kurganov-Tadmor')
#plt.plot(xe,ye,'r-',label = 'exata')
plt.xlabel( 'x' )
plt.ylabel ( 'u' )
plt.legend(loc='best')
diretorio = ""
nomefig = "problema02_b2.png"
plt.savefig(diretorio+nomefig, dpi=200)
plt.show()
