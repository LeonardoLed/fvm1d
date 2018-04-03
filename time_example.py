#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22/03/2018
@author: Leonardo Ledesma Domínguez

Example 5.1 from Malalasekera Book
----------------------------------


___________________________________________________________

 |<-------------- 1.0 m ------------->|
                                     
            u ---> 
 |------------------------------------|
                                 
\phi_0 = 1                       \phi_L = 0
___________________________________________________________



"""

import FiniteVolumeMethod as fvm
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

def analyticSol(x,t):
    return 0.5*( erfc( (x-t)/(2*(Gamma*t)**0.5) ) + np.exp(u*x/Gamma)*erfc( (x+t)/(2*(Gamma*t)**0.5) )  )

L = 2.5 # m
rho = 1.0 # kg/m^3
u = 1. # m/s
Gamma = 0.1 # kg / m.s
phi0 = 1 #
phiL = 0 #
N = 51 # Número de nodos
t_max=1.0
#cambiar esta delta a 0.0002 cuando sean 351 nodos
dt=0.002
n_tiempos=t_max/dt+1
tiempos=np.arange(0,t_max+dt,dt)
n_tiempos=len(tiempos)


#
# Creamos la malla y obtenemos datos importantes
#
malla = fvm.Mesh(nodes = N, length = L)
nx    = malla.nodes()     # Número de nodos
nvx   = malla.volumes()   # Número de volúmenes
delta = malla.delta()     # Tamaño de los volúmenes
#
# Imprimimos los datos del problema (nicely)
#
fvm.printData(Longitud = L,
              Densidad = rho,
              Velocidad = u,
              Coef_Diff = Gamma,
              Prop_0 = phi0,
              Prop_L = phiL,
              Nodos = nx, 
              Volúmenes = nvx,
              Delta = delta)
#
# Se aloja memoria para los coeficientes
#
coef = fvm.Coefficients()
coef.alloc(nvx)
#fvm.Coefficients.alloc(nvx)
#
#  Calculamos los coeficientes de FVM de la Difusión
#
dif = fvm.Diffusion1D(nvx, Gamma = Gamma, dx = delta)
dif.calcCoef()
#print('aW = {}'.format(dif.aW()), 
#      'aE = {}'.format(dif.aE()), 
#      'Su = {}'.format(dif.Su()), 
#      'aP = {}'.format(dif.aP()), sep='\n')
#print('.'+'-'*70+'.')
#
#  Calculamos los coeficientes de FVM de la Advección
#
adv = fvm.Advection1D(nvx, rho = rho, dx = delta)
adv.setU(u)
adv.calcCoef()
#print('aW = {}'.format(dif.aW()), 
#      'aE = {}'.format(dif.aE()), 
#      'Su = {}'.format(dif.Su()), 
#      'aP = {}'.format(dif.aP()), sep='\n')
#print('u = {}'.format(adv.u()))
#print('.'+'-'*70+'.')
#
# Se construye el arreglo donde se guardará la solución
#
phi = np.zeros(nvx)  # El arreglo contiene ceros
phi[0]  = phi0       # Condición de frontera izquierda
phi[-1] = phiL       # Condición de frontera derecha



# Se aplican las condiciones de frontera
#
coef.bcDirichlet('LEFT_WALL', phi0)   # Se actualizan los coeficientes
coef.bcDirichlet('RIGHT_WALL', phiL) # de acuerdo a las cond. de frontera
print('aW = {}'.format(dif.aW()), 
      'aE = {}'.format(dif.aE()), 
      'Su = {}'.format(dif.Su()), 
      'aP = {}'.format(dif.aP()), sep='\n')
print('u = {}'.format(adv.u()))
print('.'+'-'*70+'.')

# Se construye el sistema lineal de ecuaciones a partir de los coef. de FVM
#
Su = coef.Su()  # Vector del lado derecho

###--------------Definimos los indices-------------
aP=dif.aP()
aE=dif.aE()
aW=dif.aW()
Su=dif.Su()
aE[-2]=0
aW[1]=0
#------------------------------------------------
#  INICIA FOR 
#--------------------------------------------
k=1
x1 = np.linspace(0,L,350)
x = malla.createMesh()
fig = plt.figure ()
for i in range(1,n_tiempos):
        
    t=tiempos[i]
    phi_a = analyticSol(x1,t)
      
    nphi=np.copy(phi)
    for l in range(1,len(phi)-1):
        nphi[l]= (dt/(rho*delta)) * ((-aP[l]+rho*delta/dt)*phi[l]+aE[l]*phi[l+1]+aW[l]*phi[l-1]+Su[l])
    phi = nphi


    if i == int(k*(n_tiempos-1)/10):

        #print("entro")
        plt.plot(x1,phi_a, '-', label = 'Sol. analítica %.2f' %t) 
        plt.plot(x,phi,'--o', label = 'Sol. numérica')
        plt.title('Solución de $\partial(p u \phi)/\partial x= \partial (\Gamma \partial\phi/\partial x)/\partial x$ con FVM')
        plt.xlabel('$x$ [m]')
        plt.ylabel('$\phi$ [...]')
        plt.grid()
        plt.legend()

        #plt.savefig('example04.pdf')

        k=k+1



#
# Se construye un vector de coordenadas del dominio
#
plt.show()
#
# Calculamos la solución exacta y el error
#

# Calculamos la solución exacta en una malla más fina para graficar
#


#
#  Se grafica la solución
#

