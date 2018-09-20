# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 10:33:56 2018

@author: Emerson, Lucas & Emanuel
"""

import matplotlib.pyplot as plt
import numpy as np
import math as m

def exactSol2(x,t):#Solução exata para o exemplo 2

    return np.sin(x - t)

def initialCondition(Nx,h,x,ex):

    u = np.zeros([Nx])
    for i in range(0,Nx):
        if(ex == 2):
            u[i] = m.sin(x[i])
        elif(ex == 3):
            u[i] = 0.5+ m.sin(x[i])
        elif(ex == 4):
            if(x[i] < 0):
                u[i] = 2.0
            elif(x[i] > 0):
                u[i] = -2.0
            else:
                u[i] = 0
    return u

def minMod(u,Nx,h):

    theta = 2.0
    a = np.empty(Nx)
    b = np.empty(Nx) #é o c do iury
    c = np.empty(Nx) #é o b do iury
    ux = np.empty(Nx)

    if(tipoCond == 0):
        a[0] = (u[0] - u[Nx-1])/h
        a[Nx-1] = (u[Nx-1] - u[Nx-2])/h
        b[0] = (u[1] - u[Nx-1])/(2*h)
        b[Nx-1] = (u[0] - u[Nx-2])/(2*h)
        c[0] = (u[1] - u[0])/h
        c[Nx-1] = (u[0] - u[Nx-1])/h
    elif(tipoCond == 1):
        a[0] = (u[0] - u[0])/h
        a[Nx-1] = (u[Nx-1] - u[Nx-2])/h
        b[0] = (u[1] - u[0])/(2*h)
        b[Nx-1] = (u[Nx-1] - u[Nx-2])/(2*h)
        c[0] = (u[1] - u[0])/h
        c[Nx-1] = (u[Nx-1] - u[Nx-1])/h

    for i in range(1,Nx-1):
        a[i] = (u[i] - u[i-1])/h
        b[i] = (u[i+1] - u[i-1])/(2*h)
        c[i] = (u[i+1] - u[i])/h

    for i in range(Nx):
        ux[i] = 0.5*(np.sign(a[i])+np.sign(c[i]))*min(theta*abs(a[i]),abs(b[i]),theta*abs(c[i]))
        #ux[i] = 0.5*(np.sign(a[i])+np.sign(c[i]))*min(theta*abs(a[i]),theta*abs(c[i]),0.5*abs(a[i] + c[i])) #do iury

    return ux

def f(u,ex):
    if(ex == 2):
        return u
    elif(ex == 3):
        return 0.5*u**2
    elif(ex == 4):
        return 0.25*(u**2 - 1)*(u**2 - 4)

def df(u,ex):
    if(ex == 2):
        return 1.0
    elif(ex == 3):
        return abs(u)
    elif(ex == 4):
        return u**3 - (5.0/2.0)*u



"""
# Método Kurganov-Tadmor (KT)
"""
def KT(u,Nx,h):

    ux = minMod(u,Nx,h)

    H1 = np.empty(Nx)       #Hj + 1/2
    H2 = np.empty(Nx)       #Hj - 1/2
    a = np.empty(Nx)

    """
    u_minus1 = np.empty(Nx) #Hj + 1/2
    u_plus1 = np.empty(Nx)  #Hj + 1/2
    u_minus2 = np.empty(Nx) #Hj - 1/2
    u_plus2 = np.empty(Nx)  #Hj - 1/2
    a1 = np.empty(Nx)       #Hj + 1/2
    a2 = np.empty(Nx)       #Hj - 1/2

    for i in range(1,Nx-1):
        u_plus1[i]  = u[i+1] - 0.5*h*ux[i+1]
        u_minus1[i] = u[i] + 0.5*h*ux[i]
        u_plus2[i]  = u[i-1] + 0.5*h*ux[i-1]
        u_minus2[i] = u[i] - 0.5*h*ux[i]

    #Tratamento das extremidades: Circular
    u_plus1[0] = u[1] - 0.5*h*ux[1]
    u_plus2[0] = u[Nx-1] + 0.5*h*ux[Nx-1]
    u_minus1[0] = u[0] + 0.5*h*ux[0]
    u_minus2[0] = u[0] - 0.5*h*ux[0]
    u_plus1[Nx-1] = u[0] - 0.5*h*ux[0]
    u_plus2[Nx-1] = u[Nx-2] + 0.5*h*ux[Nx-2]
    u_minus1[Nx-1] = u[Nx-1] + 0.5*h*ux[Nx-1]
    u_minus2[Nx-1] = u[Nx-1] - 0.5*h*ux[Nx-1]

    for i in range(Nx):
        a1[i] = max( abs(df(u_plus1[i])) , abs(df(u_minus1[i])) )
        a2[i] = max( abs(df(u_plus2[i])) , abs(df(u_minus2[i])) )

    for i in range(0,Nx-1):
        H1[i] = ( f(u_plus1[i+1]) + f(u_minus1[i+1]))/2.0 - a1[i+1]*(u_plus1[i+1]-u_minus1[i+1])/2.0

    #Tratamento circular
    H1[Nx-1] = (f(u_plus1[0])+f(u_minus1[0]))/2.0 - a1[0]*(u_plus1[0]-u_minus1[0])/2.0

    for i in range(Nx):
        H2[i] = (f(u_plus2[i])+f(u_minus2[i]))/2.0 - a2[i]*(u_plus2[i]-u_minus2[i])/2.0
    """


    #Implementação abaixo do Iury, ainda nao revisada!
    u_minus = np.empty(Nx)
    u_plus = np.empty(Nx)

    # u+1/2
    if(tipoCond == 0):
        u_minus[0] = u[Nx-1] + h*ux[Nx-1]/2.0
    elif(tipoCond == 1):
        u_minus[0] = u[0] + h*ux[0]/2.0

    for i in range(1,Nx):
        u_minus[i] = u[i-1] + h*ux[i-1]/2.0

    # u-1/2
    for i in range(Nx):
        u_plus[i] = u[i] - h*ux[i]/2.0

    for i in range(Nx):
        a[i] = max( abs(df(u_plus[i],ex)) , abs(df(u_minus[i],ex)) )

    for i in range(0,Nx-1):
        H1[i] = ( f(u_plus[i+1],ex) + f(u_minus[i+1],ex))/2.0 - a[i+1]*(u_plus[i+1]-u_minus[i+1])/2.0


    if(tipoCond == 0):#Tratamento circular
        H1[Nx-1] = (f(u_plus[0],ex)+f(u_minus[0],ex))/2.0 - a[0]*(u_plus[0]-u_minus[0])/2.0
    elif(tipoCond == 1):
        H1[Nx-1] = H1[Nx-2]

    for i in range(Nx):
        H2[i] = (f(u_plus[i],ex)+f(u_minus[i],ex))/2.0 - a[i]*(u_plus[i]-u_minus[i])/2.0

    return (H2 - H1)/h

"""
# Método Lax-Friedrichs
"""
def laxFriedrichs(uu,Nx,Nt,h,dt):

    global K

    #Será criado 2 pontos extras que repetem o valor do contorno para tratar
    #o problema de índices nos extremos (TRATAMENTO CIRCULAR)
    u = np.empty([Nx+2,Nt+1])

    u[1:Nx+1,:] = uu[:,:]
    if(tipoCond == 0):
        u[0,0] = uu[-1,0]
        u[-1,0] = uu[0,0]
    elif(tipoCond == 1):#TRATAMENTO DUPLICAR EXTREMOS
        u[0,0] = uu[0,0]
        u[-1,0] = uu[-1,0]

    for t in range(0,Nt):
        for i in range(1,Nx+1):
            if(ex == 2):
                u[i,t+1] = 0.5*(u[i+1,t]+u[i-1,t]) - (K*dt/(2.0*h))*(u[i+1,t] - u[i-1,t])
            elif(ex == 3):
                u[i,t+1] = 0.5*(u[i+1,t]+u[i-1,t]) - (K*dt/(2.0*h))*(0.5*u[i+1,t]**2 -0.5*u[i-1,t]**2) #Exemplo 3
            elif(ex == 4):
                u[i,t+1] = 0.5*(u[i+1,t]+u[i-1,t]) - (K*dt/(2.0*h)) * 0.25 * ((u[i+1,t]**2 - 1) * (u[i+1,t]**2 - 4) - (u[i-1,t]**2 - 1) * (u[i-1,t]**2 - 4) ) #Exemplo 4
        #Copia os valores do contorno para os pontos extras criados
        if(tipoCond == 0):
            u[0,t+1] = u[Nx,t+1]
            u[Nx+1,t+1] = u[1,t+1]
        elif(tipoCond == 1):#TRATAMENTO DUPLICAR EXTREMOS
            u[0,t+1] = u[1,t+1]
            u[Nx+1,t+1] = u[Nx,t+1]

    return u[1:Nx+1,:]


ex = 4 # Exemplo 2 ou 3 ou 4
tipoCond = 1 # 0 -> Circular | 1 -> Duplicado
#---------------------------------------- Main -------------------------------#
if __name__ =="__main__":

    A = -1                                     #Limite inferior para x
    B = 1                              #Limite superior para x
    Nx2 = 160                                 #Quantidade de pontos (qtd elementos=Nx-1)
    h2 = (B - A)/np.float(Nx2-1)              #Discretização no espaço
    T = 1.2                                     #Tempo final
    xap2 = np.linspace(A,B,Nx2)               #Para plot das aproximações
    K = 1.0                                   #Coeficiente convectivo
    xe = np.linspace(A,B,1000)                #Para plot da solução exata

    u_ini2 = initialCondition(Nx2,h2,xap2,ex) #Condicao inicial


    dt2 = 0.25*h2/max(u_ini2)                 #Discretização no tempo KT aproximada
    dt3 = 1.0*h2/max(u_ini2)                 #Discretização no tempo LxF aproximada

    Nt2 = np.int(T/dt2)                       #Quantidade de iterações no tempo
    Nt3 = np.int(T/dt3)                       #Quantidade de iterações no tempo
    tj = 1.2                                  #Instante de tempo desejado
    t2 = int(tj/dt2)                          #Índice correspondente (Caso u = Matriz)
    t3 = int(tj/dt3)
    #ue = exactSol2(xe,tj)                    #Solução exata para plot
    #ue2 = exactSol2(xap,tj)                  #Solução exata para calculo do erro
    #MÉTODO KT pra exata
    if(ex == 2):
        ue = exactSol2(xe,tj)
    elif(ex == 3 or ex == 4):
        Nx = 640
        h = (B - A)/np.float(Nx-1)
        xap = np.linspace(A,B,Nx)
        u_ini = initialCondition(Nx,h,xap,ex) #Condicao inicial
        dt = 0.25*h/max(u_ini)                #Discretização no tempo KT exata
        Nt = np.int(T/dt)
        t = int(tj/dt)

        time = 0
        while(time < tj):

            #Calcula derivada de u
            du = KT(u_ini,Nx,h)

            #Realiza integração de primeira ordem para encontrar u
            uap_kt = u_ini + dt * du

            #Atribui o valor da exata para o contorno
            #uap_kt[0] = exactSol2(xap[0],time)
            #uap_kt[Nx-1] = exactSol2(xap[Nx-1],time)

            u_ini = uap_kt
            time = time + dt


    #error_kt = np.linalg.norm(np.abs(ue2-uap_kt),ord = np.inf)
    #print ('Erro KT: %f') %error_kt

    #Método KT pra aprox
    time = 0
    while(time < tj):

        #Calcula derivada de u
        du = KT(u_ini2,Nx2,h2)

        #Realiza integração de primeira ordem para encontrar u
        uap_kt2 = u_ini2 + dt2 * du

        #Atribui o valor da exata para o contorno
        #uap_kt[0] = exactSol2(xap[0],time)
        #uap_kt[Nx-1] = exactSol2(xap[Nx-1],time)

        u_ini2 = uap_kt2
        time = time + dt2

    #MÉTODO LXF
    u = np.empty([Nx2,Nt3+1])
    u[:,0] = initialCondition(Nx2,h2,xap2,ex)
    print(Nt3)
    u_lxF = laxFriedrichs(u,Nx2,Nt3,h2,dt3)
    #error_lxF = np.linalg.norm(np.abs(ue2-u_lxF[:,t]),ord = np.inf)
    #print error_lxF

    #Plota o grafico
    plt.title('Exata x Aproximada ('+ str(Nx2) +  ' pontos ,KT -> dt = ' + str(dt2)+ ', LxF -> dt = ' + str(dt3)+')')
    #plt.xlim (-0.5 ,7)
    plt.ylim (-2.2 ,2.2)
    plt.grid()

    if(ex == 2):
        plt.plot(xe,ue,'r-',label='Exata')
    elif(ex == 3 or ex == 4):
        plt.plot(xap,uap_kt,'r-',label = 'Exata')

    plt.plot(xap2,u_lxF[:,t3],'k.',label = 'Lax-Friedrichs')
    plt.plot(xap2,uap_kt2,'b.',label = 'Kurganov-Tadmor')

    plt.xlabel( 'x' )
    plt.ylabel ( 'u' )
    plt.legend(loc='best')
    diretorio = ""
    nomefig = str(Nx2)+'-lxfExemplo4.png'
    plt.savefig(diretorio+nomefig, dpi=200)
    plt.show()
