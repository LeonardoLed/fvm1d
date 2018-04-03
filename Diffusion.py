#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified on 22/03/2018
by Leonardo Ledesma Domínguez
@author:Luis Miguel de la Cruz Salas
"""


from Coefficients import Coefficients

class Diffusion1D(Coefficients):
    
    def __init__(self, nvx = None, Gamma = None, dx = None):
        super().__init__(nvx, dx)
        self.__nvx = nvx
        self.__Gamma = Gamma
        self.__dx = dx

    def __del__(self):
        del(self.__Gamma)
        del(self.__dx)
    
    def calcCoef(self):
        aE = self.aE()
        aW = self.aW()
        aP = self.aP()
        
        aE += self.__Gamma / self.__dx
        aW += self.__Gamma / self.__dx
        aP += aE + aW

        return aE,aW,aP
 
#        for i in range(self.__nvx):
#            aE[i] += self.__Gamma / self.__dx
#            aW[i] += self.__Gamma / self.__dx
#            aP[i] += aE[i] + aW[i]

if __name__ == '__main__':
    
    df1 = Diffusion1D(5, 5, 1)
    df1.alloc(5)
    df1.calcCoef()
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    df1.setSu(100)

    print('-' * 20)  
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  

    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  
