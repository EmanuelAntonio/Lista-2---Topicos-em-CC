# -*- coding: utf-8 -*-
import numpy as np
import sympy as sy

flowFunction = 0
difFlowFunction = 0


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

def readVelocityField(Nx,Ny):

    f = open('vel2.dat')
    velocityListU = []
    velocityListV = []
    velocityMatrixU = []
    velocityMatrixV = []

    contX = 0
    contY = 0
    for i in f:
        vel = i.split(',')
        l = []
        l.append(float(vel[0]))
        l.append(float(vel[1]))
        velocityListU.append(l)
        l = []
        l.append(float(vel[2]))
        l.append(float(vel[3]))
        velocityListV.append(l)
        contX += 1
        if contX == Nx:
            contX = contX % Nx
            contY += 1
            velocityMatrixU.append(velocityListU)
            velocityMatrixV.append(velocityListV)
            velocityListU = []
            velocityListV = []


    return np.array(velocityMatrixU), np.array(velocityMatrixV)




def initialCondition(Nx,Ny,h,x,y):

    sw0 = 0.21
    u = np.full((Ny,Nx),sw0)

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
            a[Ny-1,j] = (u[Ny-1,j] - u[Ny-2,j])/dy
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

def f(sw):
    global flowFunction
    """u = sy.Symbol('u') #Usa sympy para criar a funcao

    srw = 0.2
    sro = 0.15
    krw = ((u - srw)/(1 - srw))**2
    kro = (1 - u/(1-sro))**2
    mio = 10.0
    miw = 1.0

    #Define a funcao de fluxo
    F = (krw * (mio + miw))/((krw + kro) * miw)

    #Usa lambdify para melhorar desempenho
    Fl = sy.lambdify([u],F,"numpy")

    return Fl(sw) """
    return flowFunction(sw)

def df(sw):
    global difFlowFunction
    """u = sy.Symbol('u') #Usa sympy para criar a funcao

    srw = 0.2
    sro = 0.15
    krw = ((u - srw)/(1 - srw))**2
    kro = (1 - u/(1-sro))**2
    mio = 10.0
    miw = 1.0

    #Define a funcao de fluxo
    F = (krw * (mio + miw))/((krw + kro) * miw)
    dF = sy.diff(F,u)

    #Usa lambdify para melhorar desempenho
    dFl = sy.lambdify([u],dF,"numpy")

    return dFl(sw)"""
    return difFlowFunction(sw)

def createFunctions():
    global flowFunction
    global difFlowFunction
    u = sy.Symbol('u') #Usa sympy para criar a funcao

    srw = 0.2
    sro = 0.15
    krw = ((u - srw)/(1 - srw))**2
    kro = (1 - u/(1-sro))**2
    mio = 10.0
    miw = 1.0

    #Define a funcao de fluxo
    dflowFunction = (krw * (mio + miw))/((krw + kro) * miw)
    flowFunction = sy.lambdify([u],dflowFunction,"numpy")
    difFlowFunction = sy.diff(dflowFunction,u)
    difFlowFunction = sy.lambdify([u],difFlowFunction,"numpy")

def KT2D(u,Nx,Ny,h,uv,vv):

    ux = calculateUx(u,Ny,Nx,h)
    uy = calculateUy(u,Ny,Nx,h)
    phi = 0.2

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
    teste = []
    for i in range(Ny):
        for j in range(Nx):
            teste.append((uv[i,j,1],uv[i,j,0]))
            a1x[i,j] = max( abs((uv[i,j,1]/phi)*df(u_plus1x[i,j])) , abs((uv[i,j,1]/phi)*df(u_minus1x[i,j])) )
            a2x[i,j] = max( abs((uv[i,j,0]/phi)*df(u_plus2x[i,j])) , abs((uv[i,j,0]/phi)*df(u_minus2x[i,j])) )

    """
        Para Y - A
    """
    for j in range(Nx):
        for i in range(Ny):
            a1y[i,j] = max( abs((vv[i,j,1]/phi)*df(u_plus1y[i,j])) , abs((vv[i,j,1]/phi)*df(u_minus1y[i,j])) )
            a2y[i,j] = max( abs((vv[i,j,0]/phi)*df(u_plus2y[i,j])) , abs((vv[i,j,0]/phi)*df(u_minus2y[i,j])) )


    """
        Para X - H
    """
    for i in range(Ny):
        for j in range(Nx):
            H1x[i,j] = ( uv[i,j,1] * f(u_plus1x[i,j]) + uv[i,j,1] * f(u_minus1x[i,j]))/(2.0*phi) - a1x[i,j]*(u_plus1x[i,j]-u_minus1x[i,j])/2.0

    for i in range(Ny):
        for j in range(Nx):
            H2x[i,j] = (uv[i,j,0] * f(u_plus2x[i,j]) + uv[i,j,0] * f(u_minus2x[i,j]))/(2.0*phi) - a2x[i,j]*(u_plus2x[i,j] - u_minus2x[i,j])/2.0

    duX = (H2x - H1x)/h


    """
        Para Y - H
    """
    for j in range(Nx):
        for i in range(Ny):
            H1y[i,j] = ( vv[i,j,1] * f(u_plus1y[i,j]) + vv[i,j,1] * f(u_minus1y[i,j]))/(2.0*phi) - a1y[i,j]*(u_plus1y[i,j]-u_minus1y[i,j])/2.0

    for j in range(Nx):
        for i in range(Ny):
            H2y[i,j] = ( vv[i,j,0] * f(u_plus2y[i,j]) + vv[i,j,0] * f(u_minus2y[i,j]))/(2.0*phi) - a2y[i,j]*(u_plus2y[i,j] - u_minus2y[i,j])/2.0

    duY = (H2y - H1y)/h

    return duX + duY,teste

def boundaryCondition(u_ini):

    u_ini[:,0] = 1.0
    return u_ini




ex = 11 # Exemplo 2 ou 3 ou 4
tipoCond = 1 # 0 -> Circular | 1 -> Duplicado
#---------------------------------------- Main -------------------------------#
if __name__ =="__main__":

    A = 0                                    #Limite inferior para x
    B = 1.0                             #Limite superior para x
    Nx = 64                                 #Quantidade de pontos (qtd elementos=Nx-1)
    Ny = 64                               #Quantidade de pontos (qtd elementos=Nx-1)
    h = (B - A)/np.float(Nx-1)              #Discretização no espaço
    T = 0.02                                    #Tempo final
    xap = np.linspace(A,B,Nx)               #Para plot das aproximações
    yap = np.linspace(A,B,Nx)

    u_ini = initialCondition(Nx,Ny,h,xap,yap)

    dt = 0.0005*h/np.max(u_ini)                 #Discretização no tempo KT aproximada

    Nt = np.int(T/dt)                       #Quantidade de iterações no tempo
    tj = T                                  #Instante de tempo desejado

    uv,vv = readVelocityField(Nx,Ny)
    createFunctions()

    #error_kt = np.linalg.norm(np.abs(ue2-uap_kt),ord = np.inf)


    """
        Método KT para aprox
    """

    time = 0
    cont = 0

    while(time < tj):

        print ("Progress: %.2f \n" % (100.0*cont/Nt))
        u_ini = boundaryCondition(u_ini)
        #Calcula derivada de u
        dS,teste = KT2D(u_ini,Nx,Ny,h,uv,vv)

        #Salva a solução em cada passo no formato VTK -------------------------
        salvarVTK(xap,yap,u_ini,"output3/ap_"+str(cont)+".vtk")

        #Realiza integração de primeira ordem para encontrar u
        uap_kt = u_ini + dt * dS

        u_ini = uap_kt

        time = time + dt
        cont = cont + 1

    #Salva a solução no passo final no formato VTK ----------------------------
    salvarVTK(xap,yap,uap_kt,"output3/ap_"+str(cont)+".vtk")
