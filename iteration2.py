# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 13:01:07 2022
damit dz kleiner gewählt werden kann, muss dt kleiner gewählt werden
"""
import assemblymethod as am
import layermethod as lm
import numpy as np
import matplotlib.pyplot as plt

"""
create 2 layers : wie kann ich das hier machen
"""

#get all vectors
ass = am.assembly(2)

dz = ass.get_dz()
k = ass.get_k()
cv= ass.get_cv()
dzsum = ass.get_dzsum()

#time discretization
dt   = 0.1
TIME = 100
t = np.arange(0,TIME+dt,dt)
cols = len(t)

# IC, create B and A 
rows = len(dz)
A = np.zeros((rows,))
B = np.zeros((rows,))

A[:] = 1

#BC,
A[0] = 0
A[-1] = 0

"""
ITERATION FORWARD-DIFFERENCE METHOD !factor must be <0.5; 
k0:k in layer
k1:k in layer below
#relevant slices
k0s = k[1:rows-1]
k1s = k[2:rows]
cv0s = cv[1:rows-1]
cv1s = cv[2:rows]
"""
# factor vectors must be the same length as the dzsum vector

factor12 = np.zeros((rows,))
factor12u= np.zeros((rows,))
factor12l= np.zeros((rows,))

#dzn[:] = np.array(dz)
factor12[1:rows-1] = np.array(((1+k[2:rows]/k[1:rows-1])/(1+cv[1:rows-1]*k[2:rows]/cv[2:rows]/k[1:rows-1]))* cv[1:rows-1] * dt/dz[1:rows-1]**2)
factor12u[1:rows-1]= np.array(2*k[2:rows]/(k[2:rows]+k[1:rows-1]))
factor12l[1:rows-1]= np.array(2*k[1:rows-1]/(k[2:rows]+k[1:rows-1]))


#iteration zeitvektoren
for j in range(0, cols):
    #plot von A, direkt GGÜ dzsum (dz ändert sich so immer genau richtig)
    plt.plot(A[:],-dzsum)
    #Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
    B[1:rows-1] = factor12[1:rows-1]*(factor12u[1:rows-1]*A[:rows-2] - 2*A[1:rows-1] + factor12l[1:rows-1]*A[2:rows]) + A[1:rows-1]  
    A = B.copy()
    
plt.show()