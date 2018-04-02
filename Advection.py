#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22/03/2018
@author: Leonardo Ledesma Domínguez
"""


import numpy as np
from Coefficients import Coefficients


class Advection1D(Coefficients):
    
    def __init__(self, nvx = None, rho = None, dx = None):
        super().__init__(nvx)
        self.__nvx = nvx
        self.__rho = rho
        self.__dx = dx
        self.__u = np.zeros(nvx-1)

    def __del__(self):
        del(self.__nvx)
        del(self.__rho)
        del(self.__dx)
        del(self.__u)

    def setU(self, u):
        if type(u) == float:
            self.__u.fill(u)
        else:
            self.__u = u

    def u(self):
        return self.__u
    
    def calcCoef(self, method= 'CD'):
        aE = self.aE()
        aW = self.aW()
        aP = self.aP()
        u = self.__u
        rho = self.__rho

        # para métodos de segundo y tercer orden
        aEE = self.aEE()
        aWW = self.aWW()

        print ("vector de U")
        print (u)
        for i in range(1,self.__nvx-1):

            # Diferencias Centrales
            if method == 'CD':
                CE = - rho * u[i] * 0.5
                CW =   rho * u[i-1] * 0.5
                aE[i] += CE
                aW[i] += CW
                aP[i] += CE + CW + rho * (u[i] - u[i - 1])
            # Upwind de  primer orden
            elif method == 'UD':
                CE = max((-u[i],0))
                CW = max((u[i-1],0))
                aE[i] += CE
                aW[i] += CW
                aP[i] += CE + CW + rho * (u[i] - u[i - 1])
            # Upwind de segundo orden
            elif method == 'LUD':
                CE = ((3/2) * max(-rho*u[i],0)) + ((1/2) * max(-rho*u[i-1],0))
                CW = ((3/2) * max(rho*u[i-1],0)) + ((1/2) * max(rho*u[i],0))
                CEE = -(1/2) * max(-rho*u[i],0)
                CWW = -(1/2) * max(rho*u[i-1],0)
                aE[i] += CE
                aW[i] += CW
                aWW[i] += CWW
                aEE[i] += CEE
                aP[i] += CE + CW + CEE + CWW +  rho * (u[i] - u[i - 1])

            # QUICK
            elif method == 'QK':
                if (rho * u[i-1]) > 0:
                    alphaw = 1
                if (rho * u[i]) > 0 :
                    alphae = 1
                if (rho * u[i-1]) < 0 :
                    alphaw = 0
                if (rho * u[i]) < 0:
                    alphae = 0

                CE = -((3 / 8) * alphae * rho * u[i]) - ((6 / 8) * (1-alphae) * rho * u[i]) - ((1 / 8) * (1-alphaw) * rho * u[i-1])
                CW = ((6 / 8) * alphaw * rho * u[i-1]) + ((1 / 8) * alphae * rho * u[i]) + ((3 / 8) * (1-alphaw) * rho * u[i-1])
                CEE = (1 / 8) * (1-alphae) * rho * u[i]
                CWW = -(1 / 8) * alphaw * rho * u[i-1]
                aE[i] += CE
                aW[i] += CW
                aWW[i] += CWW
                aEE[i] += CEE
                aP[i] += CE + CW + CEE + CWW + rho * (u[i] - u[i - 1])



if __name__ == '__main__':
    
    nx = 5
    u = np.sin(np.linspace(0,1,nx))
#    u = np.ones(nx)
    print('-' * 20)  
    print(u)
    print('-' * 20)  

    af1 = Advection1D(6, 1, 1)
    af1.alloc(6)
    af1.setU(u)
    print(af1.u())
    print('-' * 20)  

    af1.calcCoef('UD')
    print(af1.aP(), af1.aE(), af1.aW(), af1.aWW(), af1.aEE(), af1.Su(), sep = '\n')
    print('-' * 20)  

    af1.bcDirichlet('LEFT_WALL', 2)
    af1.bcDirichlet('RIGHT_WALL', 1)
    print(af1.aP(), af1.aE(), af1.aW(), af1.aWW(), af1.aEE(),af1.Su(), sep = '\n')
    print('-' * 20)  



