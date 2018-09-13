# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import math as m

def exata(x):
    y = np.zeros(len(x))
    for i in range(len(x)):
        if(x[i] <= 0.55):
            y[i] = 1.0
    return y

def laxFriedrichs(uu,Nx,Nt,h,dt):
    global K
    u = np.zeros([Nx+1,Nt+1])
    u[:,:] = uu[:,:]
    for t in range(0,Nt):
        for i in range(1,Nx):
            u[i,t+1] = 0.5*(u[i+1,t]+u[i-1,t]) - (K*dt/(2*h))*(u[i+1,t]-u[i-1,t])
    return u

def NessyahuTadmor(uu,Nx,Nt,h,dt):
    global K
    u = np.zeros([Nx+1,Nt+1])
    u[:,:] = uu[:,:]
    for t in range(0,Nt):
        for i in range(1,Nx):
            u[(i+1)/2,t+1] = 0.5*(u[i,t]+u[i+1,t]) - (K*dt/h)*(u[i+1,(t+1)/2]-u[i-1,(t+1)/2])

def condicaoInicial(Nx,Nt,h):
    u = np.zeros([Nx+1,Nt+1])
    for i in range(0,Nx+1):
        u[i,0] = m.sin(i*h)
    return u


def condicoesContorno(u,Nx,Nt):
    for t in range(0, Nt+1):
        u[0,t] = 0.0
        u[Nx,t] = m.sin(Nx*h)
    return u

#Entrada
A = 0 #Limite inferior para x
B = m.pi #Limite superior para x
Nx = 100 #Quantidade de elementos
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

tj = 1 #Instante de tempo desejado
#t = int(tj/dt) #Índice correspondente
t = 0
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
