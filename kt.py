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

def minMod(a,b,c):
    
    if(a > 0 and b > 0 and c > 0):
        return min(a,b,c)
    elif(a < 0 and b < 0 and c < 0):
        return max(a,b,c)
    else:
        return 0

def calculateUx(u,Nx,h,teste):

    theta = 2.0
    a = np.empty(Nx)
    b = np.empty(Nx) #é o c do iury
    c = np.empty(Nx) #é o b do iury
    ux = np.empty(Nx)
    
    
    #LIXO de contorno q o emerson fez
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
    
    
    #CONTORNO de acordo com slide
    """
    if(tipoCond == 0):
        a[0] = (u[1] - u[0])/h
        a[Nx-1] = (u[0] - u[Nx - 1])/h
        
        b[0] = (u[1] - u[Nx - 1] + 2 * u[0])/(2*h)
        b[Nx-1] = (u[0] - u[Nx - 2] + 2 * u[Nx - 1])/(2*h)
        
        c[0] = (u[0] - u[Nx - 1])/h
        c[Nx - 1] = (u[Nx - 1] - u[Nx - 2])/h
    elif(tipoCond == 1):
        a[0] = (u[1] - u[0])/h
        a[Nx - 1] = (u[Nx - 1] - u[Nx - 1])/h
        
        b[0] = (u[1] - u[0] + 2 * u[0])/(2*h)
        b[Nx - 1] = (u[Nx - 1] - u[Nx - 2] + 2 * u[Nx - 2])/(2*h)
        
        c[0] = (u[0] - u[0])/h
        c[Nx - 1] = (u[Nx - 1] - u[Nx - 2])/h
    
    for i in range(1,Nx-1):
        a[i] = (u[i + 1] - u[i])/h
        b[i] = (- u[i + 1] - u[i - 1] + 2 * u[i])/(2*h)
        c[i] = (u[i] - u[i - 1])/h
    """
    
    for i in range(Nx):
        ux[i] = minMod(theta * a[i], b[i], theta * c[i])

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
    
    ux = calculateUx(u,Nx,h,1)
    
    H1 = np.empty(Nx)       #Hj + 1/2
    H2 = np.empty(Nx)       #Hj - 1/2
    a1 = np.empty(Nx)       #Hj + 1/2
    a2 = np.empty(Nx)       #Hj - 1/2
    
    u_minus1 = np.empty(Nx) #Hj + 1/2
    u_plus1 = np.empty(Nx)  #Hj + 1/2
    u_minus2 = np.empty(Nx) #Hj - 1/2
    u_plus2 = np.empty(Nx)  #Hj - 1/2


    #U + 1/2
    if(tipoCond == 0):
        u_plus1[Nx-1] = u[0] - h*ux[0]/2.0
    elif(tipoCond == 1):
        u_plus1[Nx-1] = u[Nx-1] - h*ux[Nx-1]/2.0
    
    for i in range(0,Nx-1):
        u_plus1[i]  = u[i+1] - h*ux[i+1]/2.0

    for i in range(0,Nx):
        u_minus1[i] = u[i] + h*ux[i]/2.0
        
    #U - 1/2
    for i in range(0,Nx):
        u_plus2[i]  = u[i] - h*ux[i]/2.0
    
    if(tipoCond == 0):
        u_minus2[0] = u[Nx-1] + h*ux[Nx-1]/2.0
    elif(tipoCond == 1):
        u_minus2[0] = u[0] + h*ux[0]/2.0

    for i in range(1,Nx):
        u_minus2[i] = u[i-1] + h*ux[i-1]/2.0

    
    for i in range(Nx):
        a1[i] = max( abs(df(u_plus1[i],ex)) , abs(df(u_minus1[i],ex)) )
        a2[i] = max( abs(df(u_plus2[i],ex)) , abs(df(u_minus2[i],ex)) )
    

    for i in range(0,Nx):
        H1[i] = ( f(u_plus1[i],ex) + f(u_minus1[i],ex))/2.0 - a1[i]*(u_plus1[i]-u_minus1[i])/2.0


    for i in range(Nx):
        H2[i] = (f(u_plus2[i],ex)+f(u_minus2[i],ex))/2.0 - a2[i]*(u_plus2[i]-u_minus2[i])/2.0

    return (H2 - H1)/h




ex = 4 # Exemplo 2 ou 3 ou 4
tipoCond = 1 # 0 -> Circular | 1 -> Duplicado
#---------------------------------------- Main -------------------------------#
if __name__ =="__main__":

    A = -1                                     #Limite inferior para x
    B = 1                             #Limite superior para x
    Nx2 = 80                                 #Quantidade de pontos (qtd elementos=Nx-1)
    h2 = (B - A)/np.float(Nx2-1)              #Discretização no espaço
    T = 1.2                                     #Tempo final
    xap2 = np.linspace(A,B,Nx2)               #Para plot das aproximações
    K = 1.0                                   #Coeficiente convectivo
    xe = np.linspace(A,B,1000)                #Para plot da solução exata

    u_ini1 = initialCondition(Nx2,h2,xap2,ex) #Condicao inicial
    dt2 = 0.25*h2/max(u_ini1)                 #Discretização no tempo KT aproximada

    Nt2 = np.int(T/dt2)                       #Quantidade de iterações no tempo
    tj = 1.2                                  #Instante de tempo desejado
    t2 = int(tj/dt2)                          #Índice correspondente (Caso u = Matriz)

    """
        Método KT para EXATA
    """
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
    

    """
        Método KT para aprox
    """
    time = 0
    while(time < tj):
        
        #Calcula derivada de u
        du1 = KT(u_ini1,Nx2,h2)

        #Realiza integração de primeira ordem para encontrar u
        uap_kt1 = u_ini1 + dt2 * du1
        
        u_ini1 = uap_kt1
        time = time + dt2
    
    
    #Plota o grafico
    plt.title('Exata x Aproximada ('+ str(Nx2) +  ' pontos ,KT -> dt = ' + str(dt2)+')')
    plt.ylim (-2.2 ,2.2)
    plt.grid()

    if(ex == 2):
        plt.plot(xe,ue,'r-',label='Exata')
    elif(ex == 3 or ex == 4):
        plt.plot(xap,uap_kt,'r-',label = 'Exata')

    plt.plot(xap2,uap_kt1,'b.',label = 'Kurganov-Tadmor')

    plt.xlabel( 'x' )
    plt.ylabel ( 'u' )
    plt.legend(loc='best')
    diretorio = ""
    nomefig = str(Nx2)+'-KT.png'
    plt.savefig(diretorio+nomefig, dpi=200)
    plt.show()
