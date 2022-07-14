from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt


k = 2.   # coeficiente del muelle
alpha = 1. # coeficiente de resistencia del aire
m = 1. # massa del objeto

def campo_vectorial(t,y):
    ''' Campo vectorial que define la dinámica de un oscilador armónico
         y = [x_1,x_2]'''
    der_x = y[1]
    der_y = -(k/m)*y[0] -(alpha/m)*y[1]
    return np.array([der_x,der_y])


def Sigma(x_2):
    return x_2

def dSigma(x_1,x_2):
    return -(k/m)*x_1 -(alpha/m)*x_2
# condición inicial
x_10 = 2.
x_20 = 0.
iniCond = np.array([x_10 , x_20])

#elegimos la sección tranversal y=0

#parametros
maxPaso = 1e6  #número máximo de pasos
maxInters = 100  #número máximo de intersecciones

numPaso = 1    # número de pasos
numInters = 0  #número de intersecciones
solution = integrate.DOP853(campo_vectorial, 0. , iniCond, t_bound = maxPaso, max_step = 1e-2, atol = 1e-12)
oldx_1 , oldx_2 = x_10 , x_20

tiempo = np.array([0.])
plt.figure()

while numInters < maxInters and numPaso < maxPaso:

    solution.step()
    x_1 , x_2 = solution.y

    tiempo = np.append(tiempo,solution.t)
    #print('soluciones = ',x_1, ' and ', x_2)

    if x_2*oldx_2 < 0 and oldx_2 >0:
        print(x_2,oldx_2)
        #break
        if np.abs(oldx_2)< 1e-10 or np.abs(x_2) < 1e-10:
            if np.abs(oldx_2)< 1e-10 :
                sigX_1,sigX_2  = oldx_1 , oldx_2
            else:
                sigX_1,sigX_2 = x_1, x_2

        else:
            # Newton method
            t_i , t_f = tiempo[-2] , tiempo[-1]
            x_ini , x_fin = oldx_2 , x_2
            x_c = x_2
            sol = np.zeros(2)
            while np.abs(x_c) > 1e-10:
                t_c = t_i + (t_f-t_i)/2.
                sol = solution.dense_output().__call__(t_c)
                x_c = sol[1]
                if x_ini*x_c <0:
                    x_fin , t_f = x_c , t_c
                    #print('t_f =', t_f)
                else:
                    x_ini , t_i = x_c , t_c
            sigX_1, sigX_2 = sol[0] , sol[1]
            print(sigX_1,sigX_2)


        numInters +=1
        print(sigX_1, ' and ', sigX_2)
        plt.plot(sigX_1,sigX_2, '.')
        plt.axis([-0.5, x_10, -1,2])
    else:
        pass
    numPaso +=1
    oldx_1 , oldx_2 = x_1 , x_2

plt.show()
