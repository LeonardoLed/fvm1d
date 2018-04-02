#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22/03/2018
@author: Leonardo Ledesma Domínguez
"""

import numpy as np

class Coefficients():
    """
    Esta clase define los arreglos principales para los coeficientes del
    metodo de Volumen Finito. Los arreglos son definidos como variables de
    clase para que sean compartidos por todos los objetos de esta clase.
    """    
    __aP = None
    __aE = None
    __aW = None
    __Su = None
    __nvx = None
    __delta = None

    """
    Estas variables de clases que vienen a continuación son producto del manejo de
    metodos de segundo y tercer orden para la solución de los términos advectivos,
    en una primera observancia.
    """

    __aEE = None
    __aWW = None

    def __init__(self, nvx = None, delta = None):
        Coefficients.__nvx = nvx
        Coefficients.__delta = delta


    @staticmethod
    def alloc(n):
        if Coefficients.__nvx:
            nvx = Coefficients.__nvx
        else:
            nvx = n
        Coefficients.__aP = np.zeros(nvx)
        Coefficients.__aE = np.zeros(nvx)
        Coefficients.__aW = np.zeros(nvx)
        Coefficients.__aWW = np.zeros (nvx)
        Coefficients.__aEE = np.zeros (nvx)
        Coefficients.__Su = np.zeros(nvx)
    
    def setVolumes(self, nvx):
        Coefficients.__nvx = nvx
        
    def setDelta(self, delta):
        Coefficients.__delta = delta
        
    def aP(self):
        return Coefficients.__aP

    def aE(self):
        return Coefficients.__aE
    
    def aW(self):
        return Coefficients.__aW

    def aEE(self):
        return Coefficients.__aEE

    def aWW(self):
        return Coefficients.__aWW
    
    def Su(self):
        return Coefficients.__Su

    """
        Definimos las condiciones de frontera para Upwind de primer orden y Diferencias
        Centradas utilizando como forma generalizada Dirichlet y Neumman
    """
    @staticmethod
    def bcDirichlet(wall, phi):
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        Su = Coefficients.__Su

        if wall == 'LEFT_WALL':
            aP[1] += aW[1]
            Su[1] += 2 * aW[1] * phi
        elif wall == 'RIGHT_WALL':
            aP[-2] += aE[-2]
            Su[-2] += 2 * aE[-2] * phi       

    @staticmethod
    def bcNeumman(wall, flux):
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        Su = Coefficients.__Su
        dx = Coefficients.__delta

        if wall == 'LEFT_WALL':
            aP[1] -= aW[1]
            Su[1] -= aW[1] * flux * dx
        elif wall == 'RIGHT_WALL':
            aP[-2] -= aE[-2]
            Su[-2] += aE[-2] * flux * dx

    """
    Definimos las condiciones de frontera para Upwind de Segundo orden
    utilizando como forma generalizada Dirichlet
    """
    @staticmethod
    def bcDirichlet_LUD(phiA, phiB, delta, gamma, dE, dW, dP, rho, u):
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        aWW = Coefficients.__aWW
        aEE = Coefficients.__aEE
        Su = Coefficients.__Su
        deltaMod = gamma / delta
        fvalue = rho * u

        """"
        # start node
        Sp = (2 * fvalue) - (2* dW[0])
        #Sp = (2 * fvalue) + (2* dW[0])
        aP[1] = dE[0] - Sp #+ (-(1/8)*fvalue)
        #aP[1] = (2 * fvalue)+ (3*dW[0])
        Su[1] = (-Sp) * phiA
        #print ("Calculo aP en frontera")
        #print (aP[1], Su[1], Sp)

        # second node
        aW[2] = ( ((3/2)* fvalue) + ((1/2) * fvalue) + dW[0] )
        Sp_sec = (1/2) * fvalue
        aP[2] = aW[2] + dE[2] - Sp_sec
        Su[2] = (-Sp_sec) * phiA

        # end node
        aW[-2] = ((3/2) * fvalue) - dW[-1]
        aWW[-2] = (-(1/2) * fvalue )
        Sp_end = (2 * dE[-1]) - fvalue
        aP[-2] = aW[-2] + aWW[-2] - Sp_end
        Su[-2] = (-Sp_end) * phiB

       """

        #if wall == 'LEFT_WALL':
        aP[1] += aW[1] + 3 * aWW[1]
        Su[1] += (2 * aW[1] + 4 * aWW[1]) * phiA
        aW[2] -= aWW[2]  # condición del segundo nodo (requerida en métodos de orden mayor a 1)
        Su[2] += (2 * aWW[2]) * phiA
        #elif wall == 'RIGHT_WALL':
        aP[-2] += aE[-2] + 3 * aEE[1]
        Su[-2] += (2 * aE[-2] + 4 * aEE[1]) * phiB
        aE[2] -= aEE[2]  # condición del penúltimo nodo (requerida en métodos de orden mayor a 1)
        Su[-3] += (2 * aEE[1]) * phiB

    """
        Definimos las condiciones de frontera para QUICK method
        utilizando como forma generalizada Dirichlet
    """
    @staticmethod
    def bcDirichlet_QK(phiA, phiB, delta, gamma, dE, dW, dP, rho, u):
        aP = Coefficients.__aP
        aE = Coefficients.__aE
        aW = Coefficients.__aW
        aWW = Coefficients.__aWW
        aEE = Coefficients.__aEE
        Su = Coefficients.__Su
        deltaMod = gamma / delta
        fvalue = rho * u

        # start node
        aE[1] = ((1/3) * deltaMod) - ((3/8) * fvalue) + dE[0]
        Sp = -(((8/3) * deltaMod)  +((2/8) * fvalue) + fvalue)
        aP[1] = aE[1] - Sp
        Su[1] = (-Sp) * phiA
        #print ("Calculo aP en frontera")
        #print (aP[1], Su[1], Sp)

        # second node
        aW[2] = ((7 / 8) * fvalue) + ((1 / 8) * fvalue) + dW[0]
        #aE[2] = (-(3/8) * 0.2) + dE[2]
        Sp_sec = (1/4) * fvalue
        aP[2] =  aW[2] + aE[2] - Sp_sec
        Su[2] = -(Sp_sec) * phiA

        #print ("Calculo aP en frontera")
        #print (aP[2], Su[2], aW[2], dW[0], aE)

        # end node
        aW[-2] = ((1/3)* deltaMod) + ((6/8)* fvalue) + dW[-1]
        Sp_end = -(((8 / 3) * deltaMod) - fvalue)
        aP[-2] = aW[-2] + aWW[-2] - Sp_end
        Su[-2] = (-Sp_end) * phiB
        #print ("Calculo aP en frontera")
        #print (aP[-2], aW[-2])


    def setSu(self, q):
        Su = Coefficients.__Su
        dx = Coefficients.__delta
        Su += q * dx
        
    def setSp(self, Sp):
        aP = Coefficients.__aP
        dx = Coefficients.__delta
        aP -= Sp * dx
        

if __name__ == '__main__':
    
    coef1 = Coefficients(6, 0.25)
    coef1.alloc(6)
    coef1.setSu(100)
    coef1.setSp(-2)
    
    print('-' * 20)  
    print(coef1.aP(), coef1.aE(), coef1.aW(), coef1.Su(), sep = '\n')
    print('-' * 20)  

    ap = coef1.aP()
    ap[2] = 25
    print(ap, coef1.aP(),sep='\n')
    print('-' * 20)  

    ae = coef1.aE()
    aw = coef1.aW()
    su = coef1.Su()
    ae.fill(5)
    aw.fill(5)
    ap.fill(10)
    #print (coef1.aP (), coef1.aE (), coef1.aW (), coef1.Su (), sep='\n')
    coef1.setSp(-2)
    coef1.bcDirichlet('LEFT_WALL', 2)
    coef1.bcNeumman('RIGHT_WALL', 1)
    print(coef1.aP(), coef1.aE(), coef1.aW(), coef1.Su(), sep = '\n')
    print('-' * 20)  

