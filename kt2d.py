import matplotlib.pyplot as plt
import numpy as np
#import math as m
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def initialCondition(Nx,Ny,h,x,y,ex):
    
    u = np.zeros([Nx,Ny])
    if(ex == 11):
        center = 1/2
        for i in range(0,Nx):
    
                if(-center - 0.4 < x[i] < -center + 0.4):
                    u[i,0] = 1.0
                elif(center - 0.4 < x[i] < center + 0.4):
                    u[i,0] = -1.0
                else:
                    u[i,0] = 0
                    
        for i in range(0,Ny):
            
                if(-center - 0.4 < y[i] < -center + 0.4):
                    u[0,i] = 1.0
                elif(center - 0.4 < y[i] < center + 0.4):
                    u[0,i] = -1.0
                else:
                    u[0,i] = 0
    return u

def minMod(a,b,c):

    if(a > 0 and b > 0 and c > 0):
        return min(a,b,c)
    elif(a < 0 and b < 0 and c < 0):
        return max(a,b,c)
    else:
        return 0

def calculateUx(u,Nx,Ny,dx):

    theta = 2.0
    a = np.empty([Nx,Ny])
    b = np.empty([Nx,Ny])
    c = np.empty([Nx,Ny])
    ux = np.empty([Nx,Ny])

    #LIXO de contorno q o emerson fez
    for k in range(Ny):
        if(tipoCond == 0):
            a[0,k] = (u[0,k] - u[Nx-1,k])/dx
            a[Nx-1,k] = (u[Nx-1,k] - u[Nx-2,k])/dx
            b[0,k] = (u[1,k] - u[Nx-1,k])/(2*dx)
            b[Nx-1,k] = (u[0,k] - u[Nx-2,k])/(2*dx)
            c[0,k] = (u[1,k] - u[0,k])/dx
            c[Nx-1,k] = (u[0,k] - u[Nx-1,k])/dx
        elif(tipoCond == 1):
            a[0,k] = (u[0,k] - u[0,k])/dx
            a[Nx-1,k] = (u[Nx-1,k] - u[Nx-2,k])/dx
            b[0,k] = (u[1,k] - u[0,k])/(2*dx)
            b[Nx-1,k] = (u[Nx-1,k] - u[Nx-2,k])/(2*dx)
            c[0,k] = (u[1,k] - u[0,k])/dx
            c[Nx-1,k] = (u[Nx-1,k] - u[Nx-1,k])/dx

    for k in range(Ny):
        for i in range(1,Nx-1):
            a[i,k] = (u[i,k] - u[i-1,k])/dx
            b[i,k] = (u[i+1,k] - u[i-1,k])/(2*dx)
            c[i,k] = (u[i+1,k] - u[i,k])/dx
    
    
    for k in range(Ny):
        for i in range(Nx):
            ux[i,k] = minMod(theta * a[i,k], b[i,k], theta * c[i,k])


    return ux

def calculateUy(u,Nx,Ny,dy):

    theta = 2.0
    a = np.empty([Nx,Ny])
    b = np.empty([Nx,Ny])
    c = np.empty([Nx,Ny])
    uy = np.empty([Nx,Ny])

    #LIXO de contorno q o emerson fez
    for k in range(Nx):
        if(tipoCond == 0):
            a[k,0] = (u[k,0] - u[k,Nx-1])/dy
            a[k,Nx-1] = (u[k,Nx-1] - u[k,Nx-2])/dy
            b[k,0] = (u[k,1] - u[k,Nx-1])/(2*dy)
            b[k,Nx-1] = (u[k,0] - u[k,Nx-2])/(2*dy)
            c[k,0] = (u[k,1] - u[k,0])/dy
            c[k,Nx-1] = (u[k,0] - u[k,Nx-1])/dy
        elif(tipoCond == 1):
            a[k,0] = (u[k,0] - u[k,0])/dy
            a[k,Nx-1] = (u[k,Nx-1] - u[k,Nx-2])/dy
            b[k,0] = (u[k,1] - u[k,0])/(2*dy)
            b[k,Nx-1] = (u[k,Nx-1] - u[k,Nx-2])/(2*dy)
            c[k,0]= (u[k,1] - u[k,0])/dy
            c[k,Nx-1] = (u[k,Nx-1] - u[k,Nx-1])/dy
            
    for k in range(Nx):
        for i in range(1,Nx-1):
            a[k,i] = (u[k,i] - u[k,i-1])/dy
            b[k,i] = (u[k,i+1] - u[k,i-1])/(2*dy)
            c[k,i] = (u[k,i+1] - u[k,i])/dy
    
    
    
    for k in range(Nx):
        for i in range(Ny):
            uy[k,i] = minMod(theta * a[k,i], b[k,i], theta * c[k,i])

    return uy

def f(u,ex):
    
    if(ex == 11):
        return u**2

def df(u,ex):

    if(ex == 11):
        return 2*u


"""
# Método Kurganov-Tadmor (KT)
"""

def KT2D(u,Nx,Ny,h):

    ux = calculateUx(u,Nx,Ny,h)
    uy = calculateUy(u,Nx,Ny,h)
    fi = 1
    
    uv = np.ones([Nx,Ny]) #CAMPO de velocidade em X FALTA PREENCHER
    vv = np.ones([Nx,Ny]) #CAMPO de velocidade em Y FALTA PREENCHER

    H1x = np.empty([Nx,Ny])       #Hj + 1/2
    H2x = np.empty([Nx,Ny])       #Hj - 1/2
    a1x = np.empty([Nx,Ny])       #Hj + 1/2
    a2x = np.empty([Nx,Ny])       #Hj - 1/2

    u_minus1x = np.empty([Nx,Ny]) #Hj + 1/2
    u_plus1x = np.empty([Nx,Ny])  #Hj + 1/2
    u_minus2x = np.empty([Nx,Ny]) #Hj - 1/2
    u_plus2x = np.empty([Nx,Ny])  #Hj - 1/2

    H1y = np.empty([Nx,Ny])       #Hj + 1/2
    H2y = np.empty([Nx,Ny])       #Hj - 1/2
    a1y = np.empty([Nx,Ny])       #Hj + 1/2
    a2y = np.empty([Nx,Ny])       #Hj - 1/2

    u_minus1y = np.empty([Nx,Ny]) #Hj + 1/2
    u_plus1y = np.empty([Nx,Ny])  #Hj + 1/2
    u_minus2y = np.empty([Nx,Ny]) #Hj - 1/2
    u_plus2y = np.empty([Nx,Ny])  #Hj - 1/2

    
    """
        Para X - U
    """
    #U+ + 1/2
    for k in range(Ny):
        if(tipoCond == 0):
            u_plus1x[Nx-1,k] = u[0,k] - h*ux[0,k]/2.0
        elif(tipoCond == 1):
            u_plus1x[Nx-1,k] = u[Nx-1,k] - h*ux[Nx-1,k]/2.0
    
    for k in range(Ny):
        for i in range(0,Nx-1):
            u_plus1x[i,k]  = u[i+1,k] - h*ux[i+1,k]/2.0
    
    #U- + 1/2
    for k in range(Ny):
        for i in range(0,Nx):
            u_minus1x[i,k] = u[i,k] + h*ux[i,k]/2.0

    #U+ - 1/2
    for k in range(Ny):
        for i in range(0,Nx):
            u_plus2x[i,k]  = u[i,k] - h*ux[i,k]/2.0
    
    #U- - 1/2
    for k in range(Ny):
        if(tipoCond == 0):
            u_minus2x[0,k] = u[Nx-1,k] + h*ux[Nx-1,k]/2.0
        elif(tipoCond == 1):
            u_minus2x[0,k] = u[0,k] + h*ux[0,k]/2.0
    
    for k in range(Ny):
        for i in range(1,Nx):
            u_minus2x[i,k] = u[i-1,k] + h*ux[i-1,k]/2.0

    """
        Para Y - U
    """
    #U+ + 1/2
    for k in range(Nx):
        if(tipoCond == 0):
            u_plus1y[k,Nx-1] = u[k,0] - h*uy[k,0]/2.0
        elif(tipoCond == 1):
            u_plus1y[k,Nx-1] = u[k,Nx-1] - h*uy[k,Nx-1]/2.0
    
    for k in range(Nx):
        for i in range(0,Ny-1):
            u_plus1y[k,i]  = u[k,i+1] - h*uy[k,i+1]/2.0
    
    #U- + 1/2
    for k in range(Nx):
        for i in range(0,Ny):
            u_minus1y[k,i] = u[k,i] + h*uy[k,i]/2.0

    #U+ - 1/2
    for k in range(Nx):
        for i in range(0,Ny):
            u_plus2y[k,i]  = u[k,i] - h*uy[k,i]/2.0
    
    #U- - 1/2
    for k in range(Nx):
        if(tipoCond == 0):
            u_minus2y[k,0] = u[k,Nx-1] + h*uy[k,Nx-1]/2.0
        elif(tipoCond == 1):
            u_minus2y[k,0] = u[k,0] + h*uy[k,0]/2.0
    
    for k in range(Nx):
        for i in range(1,Ny):
            u_minus2y[k,i] = u[k,i-1] + h*uy[k,i-1]/2.0
    
    
    
    """
        Para X - A
    """
    for k in range(Ny):
        for i in range(Nx):
            a1x[i,k] = max( abs((uv[i,k]/fi)*df(u_plus1x[i,k],ex)) , abs((uv[i,k]/fi)*df(u_minus1x[i,k],ex)) )
            a2x[i,k] = max( abs((uv[i,k]/fi)*df(u_plus2x[i,k],ex)) , abs((uv[i,k]/fi)*df(u_minus2x[i,k],ex)) )

    """
        Para Y - A
    """
    for k in range(Nx):
        for i in range(Ny):
            a1y[k,i] = max( abs((vv[i,k]/fi)*df(u_plus1y[k,i],ex)) , abs((vv[i,k]/fi)*df(u_minus1y[k,i],ex)) )
            a2y[k,i] = max( abs((vv[i,k]/fi)*df(u_plus2y[k,i],ex)) , abs((vv[i,k]/fi)*df(u_minus2y[k,i],ex)) )
    
    
    """
        Para X - H
    """
    for k in range(Ny):
        for i in range(0,Nx):
            H1x[i,k] = ( f(u_plus1x[i,k],ex) + f(u_minus1x[i,k],ex))/(2.0*fi) - a1x[i,k]*(u_plus1x[i,k]-u_minus1x[i,k])/2.0

    for k in range(Ny):
        for i in range(Nx):
            H2x[i,k] = (f(u_plus2x[i,k],ex)+f(u_minus2x[i,k],ex))/(2.0*fi) - a2x[i,k]*(u_plus2x[i,k] - u_minus2x[i,k])/2.0

    duX = (H2x - H1x)/h


    """
        Para Y - H
    """
    for k in range(Nx):
        for i in range(0,Ny):
            H1y[i,k] = ( f(u_plus1y[k,i],ex) + f(u_minus1y[k,i],ex))/(2.0*fi) - a1y[k,i]*(u_plus1y[k,i]-u_minus1y[k,i])/2.0

    for k in range(Nx):
        for i in range(Ny):
            H2y[i,k] = (f(u_plus2y[k,i],ex)+f(u_minus2y[k,i],ex))/(2.0*fi) - a2y[k,i]*(u_plus2y[k,i] - u_minus2y[k,i])/2.0
    
    duY = (H2y - H1y)/h
    
    
    return duX + duY




ex = 11 # Exemplo 2 ou 3 ou 4
tipoCond = 1 # 0 -> Circular | 1 -> Duplicado
#---------------------------------------- Main -------------------------------#
if __name__ =="__main__":
    
    A = -1.5                                    #Limite inferior para x
    B = 1.5                             #Limite superior para x
    Nx2 = 100                                 #Quantidade de pontos (qtd elementos=Nx-1)
    Ny2 = 100                                 #Quantidade de pontos (qtd elementos=Nx-1)
    h2 = (B - A)/np.float(Nx2-1)              #Discretização no espaço
    T = 0.5                                     #Tempo final
    xap2 = np.linspace(A,B,Nx2)               #Para plot das aproximações
    yap2 = np.linspace(A,B,Nx2)
    
    u_ini1 = initialCondition(Nx2,Ny2,h2,xap2,yap2,ex) #Condicao inicial
    dt2 = 0.25*h2/max(u_ini1[:,0])                 #Discretização no tempo KT aproximada

    Nt2 = np.int(T/dt2)                       #Quantidade de iterações no tempo
    tj = 0.5                                  #Instante de tempo desejado
    t2 = int(tj/dt2)                          #Índice correspondente (Caso u = Matriz)
    
    #error_kt = np.linalg.norm(np.abs(ue2-uap_kt),ord = np.inf)
    
    
    """
        Método KT para aprox
    """
    
    time = 0
    while(time < tj):

        #Calcula derivada de u
        du1 = KT2D(u_ini1,Nx2,Ny2,h2)

        #Realiza integração de primeira ordem para encontrar u
        uap_kt1 = u_ini1 + dt2 * du1

        u_ini1 = uap_kt1
        time = time + dt2

    
    #Plota o grafico
    print(type(uap_kt1))
    print(uap_kt1)
    mesh = Poly3DCollection(uap_kt1)
    
    #mesh.set_edgecolor("yellow")
    fig = plt.figure(figsize=(20, 20),edgecolor="g")
    ax = fig.add_subplot(111, projection='3d')
    
    ax.add_collection3d(mesh)
    
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    
    ax.set_xlim(- 2, + 2)
    ax.set_ylim(- 2, + 2)
    ax.set_zlim(- 2, + 2)
    
    plt.show()