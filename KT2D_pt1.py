# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


#Salva a solução no formato VTK -----------------------------------------------
def salvarVTK(x,y,u,diretorio):
    nx = len(x)
    ny = len(y)
    #Cabeçalho
    out = open(diretorio, 'w')
    out.write('# vtk DataFile Version 2.0\n')
    out.write('vtk output\n')
    out.write('ASCII\n')
    out.write('DATASET STRUCTURED_GRID\n')
    out.write('DIMENSIONS '+str(nx)+' '+str(ny)+' 1\n')
    out.write('POINTS '+str(nx*ny)+' float\n')
    #Coordenadas (varia primeiro em x, depois em y)
    for i in range(ny):
        for j in range(nx):
            a = u[i,j]#*10000 #Para plotar elevação
            out.write(str(x[j])+' '+str(y[i])+' '+str(a)+'\n')

    #Solução
    out.write('POINT_DATA '+str(nx*ny)+'\n')
    out.write('SCALARS scalars double\n')
    out.write('LOOKUP_TABLE default\n')
    for i in range(ny):
        for j in range(nx):
            out.write(str(u[i,j])+'\n')
    out.close()


def initialCondition(Nx,Ny,h,x,y,ex):

    u = np.empty([Ny,Nx])
    if(ex == 11):
        cx = 0.5
        cy = cx
        radius = 0.4

        for i in range(Ny):
            for j in range(Nx):
                if((x[j] - cx)**2 + (y[i] - cy)**2 <= radius**2):
                    u[i,j] = -1.0
                elif((x[j] + cx)**2 + (y[i] + cy)**2 <= radius**2):
                    u[i,j] = 1.0
                else:
                    u[i,j] = 0

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
    a = np.empty([Ny,Nx])
    b = np.empty([Ny,Nx])
    c = np.empty([Ny,Nx])
    ux = np.empty([Ny,Nx])

    for i in range(Ny):
        if(tipoCond == 0):
            a[i,0] = (u[i,0] - u[i,Nx-1])/dx
            a[i,Nx-1] = (u[i,Nx-1] - u[i,Nx-2])/dx
            b[i,0] = (u[i,1] - u[i,Nx-1])/(2*dx)
            b[i,Nx-1] = (u[i,0] - u[i,Nx-2])/(2*dx)
            c[i,0] = (u[i,1] - u[i,0])/dx
            c[i,Nx-1] = (u[i,0] - u[i,Nx-1])/dx
        elif(tipoCond == 1):
            a[i,0] = (u[i,0] - u[i,0])/dx
            a[i,Nx-1] = (u[i,Nx-1] - u[i,Nx-2])/dx
            b[i,0] = (u[i,1] - u[i,0])/(2*dx)
            b[i,Nx-1] = (u[i,Nx-1] - u[i,Nx-2])/(2*dx)
            c[i,0] = (u[i,1] - u[i,0])/dx
            c[i,Nx-1] = (u[i,Nx-1] - u[i,Nx-1])/dx

    for i in range(Ny):
        for j in range(1,Nx-1):
            a[i,j] = (u[i,j] - u[i,j-1])/dx
            b[i,j] = (u[i,j+1] - u[i,j-1])/(2*dx)
            c[i,j] = (u[i,j+1] - u[i,j])/dx


    for i in range(Ny):
        for j in range(Nx):
            ux[i,j] = minMod(theta * a[i,j], b[i,j], theta * c[i,j])


    return ux

def calculateUy(u,Nx,Ny,dy):

    theta = 2.0
    a = np.empty([Ny,Nx])
    b = np.empty([Ny,Nx])
    c = np.empty([Ny,Nx])
    uy = np.empty([Ny,Nx])

    for j in range(Nx):
        if(tipoCond == 0):
            a[0,j] = (u[0,j] - u[Ny-1,j])/dy
            a[Ny-1,j] = (u[Ny-1,j] - u[Ny-2],j)/dy
            b[0,j] = (u[1,j] - u[Ny-1,j])/(2*dy)
            b[Ny-1,j] = (u[0,j] - u[Ny-2,j])/(2*dy)
            c[0,j] = (u[1,j] - u[0,j])/dy
            c[Ny-1,j] = (u[0,j] - u[Ny-1,j])/dy
        elif(tipoCond == 1):
            a[0,j] = (u[0,j] - u[0,j])/dy
            a[Ny-1,j] = (u[Ny-1,j] - u[Ny-2,j])/dy
            b[0,j] = (u[1,j] - u[0,j])/(2*dy)
            b[Ny-1,j] = (u[Ny-1,j] - u[Ny-2,j])/(2*dy)
            c[0,j]= (u[1,j] - u[0,j])/dy
            c[Ny-1,j] = (u[Ny-1,j] - u[Ny-1,j])/dy

    for j in range(Nx):
        for i in range(1,Ny-1):
            a[i,j] = (u[i,j] - u[i-1,j])/dy
            b[i,j] = (u[i+1,j] - u[i-1,j])/(2*dy)
            c[i,j] = (u[i+1,j] - u[i,j])/dy



    for j in range(Nx):
        for i in range(Ny):
            uy[i,j] = minMod(theta * a[i,j], b[i,j], theta * c[i,j])

    return uy

def f(u,ex):

    if(ex == 11):
        return u**2

def df(u,ex):

    if(ex == 11):
        return 2*u




def KT2D(u,Nx,Ny,h):

    ux = calculateUx(u,Ny,Nx,h)
    uy = calculateUy(u,Ny,Nx,h)
    fi = 1

    uv = np.ones([Ny,Nx]) #CAMPO de velocidade em X FALTA PREENCHER
    vv = np.ones([Ny,Nx]) #CAMPO de velocidade em Y FALTA PREENCHER

    H1x = np.empty([Ny,Nx])       #Hj + 1/2
    H2x = np.empty([Ny,Nx])       #Hj - 1/2
    a1x = np.empty([Ny,Nx])       #Hj + 1/2
    a2x = np.empty([Ny,Nx])       #Hj - 1/2

    u_minus1x = np.empty([Ny,Nx]) #Hj + 1/2
    u_plus1x = np.empty([Ny,Nx])  #Hj + 1/2
    u_minus2x = np.empty([Ny,Nx]) #Hj - 1/2
    u_plus2x = np.empty([Ny,Nx])  #Hj - 1/2

    H1y = np.empty([Ny,Nx])       #Hj + 1/2
    H2y = np.empty([Ny,Nx])       #Hj - 1/2
    a1y = np.empty([Ny,Nx])       #Hj + 1/2
    a2y = np.empty([Ny,Nx])       #Hj - 1/2

    u_minus1y = np.empty([Ny,Nx]) #Hj + 1/2
    u_plus1y = np.empty([Ny,Nx])  #Hj + 1/2
    u_minus2y = np.empty([Ny,Nx]) #Hj - 1/2
    u_plus2y = np.empty([Ny,Nx])  #Hj - 1/2


    """
        Para X - U
    """
    #U+ + 1/2
    #Tratamento para o contorno
    for i in range (Ny):
        if(tipoCond == 0):
            u_plus1x[i,Nx-1] = u[i,0] - h*ux[i,0]/2.0
        elif(tipoCond == 1):
            u_plus1x[i,Nx-1] = u[i,Nx-1] - h*ux[i,Nx-1]/2.0

    #Caso geral
    for i in range(Ny):
        for j in range(Nx-1):
            u_plus1x[i,j]  = u[i,j+1] - h*ux[i,j+1]/2.0

    #U- + 1/2
    for i in range(Ny):
        for j in range(Nx):
            u_minus1x[i,j] = u[i,j] + h*ux[i,j]/2.0


    #U+ - 1/2
    for i in range(Ny):
        for j in range(Nx):
            u_plus2x[i,j] = u[i,j] - h*ux[i,j]/2.0

    #U- - 1/2
    #Tratamento para o contorno
    for i in range (Ny):
        if(tipoCond == 0):
            u_minus2x[i,0] = u[i,Nx-1] + h*ux[i,Nx-1]/2.0
        elif(tipoCond == 1):
            u_minus2x[i,0] = u[i,0] + h*ux[i,0]/2.0

    #Caso geral
    for i in range(Ny):
        for j in range(1,Nx):
            u_minus2x[i,j]  = u[i,j-1] + h*ux[i,j-1]/2.0

    """
        Para Y - U
    """
    #U+ + 1/2
    #Tratamento para o contorno
    for j in range (Nx):
        if(tipoCond == 0):
            u_plus1y[Nx-1,j] = u[0,j] - h*uy[0,j]/2.0
        elif(tipoCond == 1):
            u_plus1y[Nx-1,j] = u[Nx-1,j] - h*uy[Nx-1,j]/2.0

    #Caso geral
    for j in range(Nx):
        for i in range(Ny-1):
            u_plus1y[i,j]  = u[i+1,j] - h*uy[i+1,j]/2.0

    #U- + 1/2
    for j in range(Nx):
        for i in range(Ny):
            u_minus1y[i,j] = u[i,j] + h*uy[i,j]/2.0


    #U+ - 1/2
    for j in range(Nx):
        for i in range(Ny):
            u_plus2y[i,j] = u[i,j] - h*uy[i,j]/2.0

    #U- - 1/2
    #Tratamento para o contorno
    for j in range (Nx):
        if(tipoCond == 0):
            u_minus2y[0,j] = u[Nx-1,j] + h*uy[Nx-1,j]/2.0
        elif(tipoCond == 1):
            u_minus2y[0,j] = u[0,j] + h*uy[0,j]/2.0

    #Caso geral
    for j in range(Nx):
        for i in range(1,Ny):
            u_minus2y[i,j]  = u[i-1,j] + h*uy[i-1,j]/2.0



    """
        Para X - A
    """
    for i in range(Ny):
        for j in range(Nx):
            a1x[i,j] = max( abs((uv[i,j]/fi)*df(u_plus1x[i,j],ex)) , abs((uv[i,j]/fi)*df(u_minus1x[i,j],ex)) )
            a2x[i,j] = max( abs((uv[i,j]/fi)*df(u_plus2x[i,j],ex)) , abs((uv[i,j]/fi)*df(u_minus2x[i,j],ex)) )

    """
        Para Y - A
    """
    for j in range(Nx):
        for i in range(Ny):
            a1y[i,j] = max( abs((vv[i,j]/fi)*df(u_plus1y[i,j],ex)) , abs((vv[i,j]/fi)*df(u_minus1y[i,j],ex)) )
            a2y[i,j] = max( abs((vv[i,j]/fi)*df(u_plus2y[i,j],ex)) , abs((vv[i,j]/fi)*df(u_minus2y[i,j],ex)) )


    """
        Para X - H
    """
    for i in range(Ny):
        for j in range(Nx):
            H1x[i,j] = ( f(u_plus1x[i,j],ex) + f(u_minus1x[i,j],ex))/(2.0*fi) - a1x[i,j]*(u_plus1x[i,j]-u_minus1x[i,j])/2.0

    for i in range(Ny):
        for j in range(Nx):
            H2x[i,j] = (f(u_plus2x[i,j],ex)+f(u_minus2x[i,j],ex))/(2.0*fi) - a2x[i,j]*(u_plus2x[i,j] - u_minus2x[i,j])/2.0

    duX = (H2x - H1x)/h


    """
        Para Y - H
    """
    for j in range(Nx):
        for i in range(Ny):
            H1y[i,j] = ( f(u_plus1y[i,j],ex) + f(u_minus1y[i,j],ex))/(2.0*fi) - a1y[i,j]*(u_plus1y[i,j]-u_minus1y[i,j])/2.0

    for j in range(Nx):
        for i in range(Ny):
            H2y[i,j] = (f(u_plus2y[i,j],ex)+f(u_minus2y[i,j],ex))/(2.0*fi) - a2y[i,j]*(u_plus2y[i,j] - u_minus2y[i,j])/2.0

    duY = (H2y - H1y)/h


    return duX + duY




ex = 11 # Exemplo 2 ou 3 ou 4
tipoCond = 1 # 0 -> Circular | 1 -> Duplicado
#---------------------------------------- Main -------------------------------#
if __name__ =="__main__":

    A = -1.5                                    #Limite inferior para x
    B = 1.5                             #Limite superior para x
    Nx = 60                                 #Quantidade de pontos (qtd elementos=Nx-1)
    Ny = 60                               #Quantidade de pontos (qtd elementos=Nx-1)
    h = (B - A)/np.float(Nx-1)              #Discretização no espaço
    T = 1.0                                     #Tempo final
    xap = np.linspace(A,B,Nx)               #Para plot das aproximações
    yap = np.linspace(A,B,Nx)

    u_ini = initialCondition(Nx,Ny,h,xap,yap,ex)

    dt = 0.120*h/np.max(u_ini)                 #Discretização no tempo KT aproximada

    Nt = np.int(T/dt)                       #Quantidade de iterações no tempo
    tj = 0.5                                  #Instante de tempo desejado

    #error_kt = np.linalg.norm(np.abs(ue2-uap_kt),ord = np.inf)


    """
        Método KT para aprox
    """

    time = 0
    cont = 0
    while(time < tj):

        #Calcula derivada de u
        dS = KT2D(u_ini,Nx,Ny,h)

        #Salva a solução em cada passo no formato VTK -------------------------
        salvarVTK(xap,yap,u_ini,"output3/ap_"+str(cont)+".vtk")

        #Realiza integração de primeira ordem para encontrar u
        uap_kt = u_ini + dt * dS

        u_ini = uap_kt

        time = time + dt
        cont = cont + 1

    #Salva a solução no passo final no formato VTK ----------------------------
    salvarVTK(xap,yap,uap_kt,"output3/ap_"+str(cont)+".vtk")
