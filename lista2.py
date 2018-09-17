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
    for t in range(0,Nt):
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
def kurganovTadmor(uu,Nx,Nt,h,dt):
    u = np.zeros([Nx+1,Nt+1])
    u[:,:] = uu[:,:]

    for t in range(0,Nt):
        for j in range(1,Nx-1):
            u[j,t] = - (H(u,j+0.5,t,h)-H(u,j-0.5,t,h))/h

    return u

def H(u,j,t,h):
    return 0.5*(f(uPlus(u,j,t,h)) + f(uMinus(u,j,t,h))) - (0.5*a(u,j,t,h))*(uPlus(u,j,t,h) - uMinus(u,j,t,h))

def uPlus(u,j,t,h):
    up[j,t] = u[j + 0.5,t] - 0.5*h*uX(u,j + 0.5,t,h)
    return up

def uMinus(u,j,t,h):
    um[j,t] = u[j - 0.5] + 0.5*h*uX(u,j-0.5,t,h)
    return um

def uX(u,j,t,h):
    lista = list()
    lista.append(2*(u[j,t] - u[j-1,t])/h)
    lista.append((u[j+1,t] - u[j-1,t])/2*h)
    lista.append(2*(u[j+1,t] - u[j,t])/h)

    resultado = min(lista)

    return abs(resultado)

def a(u,j,t,h):
    return max(derivadaF(uPlus(u,j,t,h)),derivadaF(uMinus(u,j,t,h)))

def f(x):
    return x

def derivadaF(x):
    return x

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
B = m.pi #Limite superior para x
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

u_laxFri = laxFriedrichs(u,Nx,Nt,h,dt)
#u_kt = kurganovTadmor(u,Nx,Nt,h,dt)

tj = 0 #Instante de tempo desejado
t = int(tj/dt) #Índice correspondente

#Plota o grafico
#plt.title('Exata x Aproximada (400 elementos, dt = 1.5h)')
plt.title('Exata x Aproximada ('+ str(Nx) +  ' elementos , dt = ' + str(dt/h) + 'h)')
#plt.xlim (-1 ,1)
#plt.ylim (0 ,1.1)
plt.grid()

plt.plot(x,u_laxFri[:,t],'b-',label = 'Lax-Friedrichs')
#plt.plot(xe,ye,'r-',label = 'exata')
plt.xlabel( 'x' )
plt.ylabel ( 'u' )
plt.legend(loc='best')
diretorio = "/home/lucas/Downloads/"
nomefig = "problema02_b2.png"
plt.savefig(diretorio+nomefig, dpi=200)
plt.show()
