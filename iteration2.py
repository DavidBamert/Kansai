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
dt   = 0.5
TIME = 1000
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
up = slice(0,rows-2)
zero = slice(1,rows-1)
lo = slice(2, rows)

# factor vectors must be the same length as the dzsum vector

factor12 = np.zeros((rows,))
factor12u= np.zeros((rows,))
factor12l= np.zeros((rows,))

#dzn[:] = np.array(dz)
factor12[zero] = np.array( ( (1+k[lo]/k[zero]) / (1+cv[zero]*k[lo]/cv[lo]/k[zero]) ) *cv[zero]*dt/dz[zero]**2 )
factor12u[zero]= np.array( 2*k[lo] / (k[lo]+k[zero]) )
factor12l[zero]= np.array( 2*k[zero] / (k[lo]+k[zero]) )


#iteration zeitvektoren
for j in range(0, cols):
    #plot von A, direkt GGÜ dzsum (dz ändert sich so immer genau richtig)
    plt.plot(A[:],-dzsum)
    #Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
    B[zero] = factor12[zero] * ( factor12u[zero]*A[up] - 2*A[zero] + factor12l[zero]*A[lo] ) + A[zero]
    A = B.copy()
    
plt.show()