# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 13:01:07 2022
damit dz kleiner gewählt werden kann, muss dt kleiner gewählt werden
"""
import assemblymethod as am
import layermethod as lm
import numpy as np
import matplotlib.pyplot as plt

#discretizazion
#damit dz kleiner gewählt werden kann, muss dt kleiner gewählt werden
dz = 0.5
dt = 0.5

#fill Layerlist with Layers, from top to bottom
#lm.layer(hup, hlow, k, me, dz), hup und hlow von oben gemessen(+)

L = [lm.Layer(0, 10, 0.2, 0.3, dz),
    lm.Layer(10, 20, 1, 0.3, dz),
    lm.Layer(20, 30, 0.2, 0.3, dz),
    ]

#throw L into Assembly method
ass = am.Assembly(L,dt)

#get all the iteration vektors
dz = ass.get_dz()
k = ass.get_k()
cv= ass.get_cv()
dzsum = ass.get_dzsum()

#time discretization
TIME = 10000
t = np.arange(0,TIME+dt,dt)
cols = len(t)

# IC, create B and A
rows = len(dz)
A = np.zeros((rows,))
B = np.zeros((rows,))

#IC
A[:] = 0

#BCs, [upper, lower] #BC muss in A und B geändert werden. (weil nur in slice 'zero' iteriert wird)
bcs = [0, 0]

A[0], B[0] = bcs[0], bcs[0]
A[-1], B[-1] = bcs[1], bcs[1]


"""
ITERATION FORWARD-DIFFERENCE METHOD !factor must be <0.5; 
up = i-1
zero = i
lo = i+1
"""

up = slice(0,rows-2)
zero = slice(1,rows-1)
lo = slice(2, rows)

# factor vectors must be the same length as the dzsum vector

fv = np.zeros((rows,))
f1= np.zeros((rows,))
f2= np.zeros((rows,))

#dzn[:] = np.array(dz)
fv[zero] = np.array(((1 + k[zero] / k[lo]) / (1 + cv[lo] * k[zero] / cv[zero] / k[lo])) * cv[lo] * dt / dz[zero] ** 2)
f1[zero]= np.array(2 * k[zero] / (k[lo] + k[zero]))
f2[zero]= np.array(2 * k[lo] / (k[lo] + k[zero]))

# complete the first and last entry
fv[0],fv[-1] = fv[1],fv[-2]
f1[0],f1[-1] = f1[1],f1[-2]
f2[0],f2[-1] = f2[1],f2[-2]

#additional loads in time tl = np.array([[time,load], ... ]) Matrix kann beliebig erweitert werden. Eintrag [0,1] kann IC ersetzen.
tl = np.array([
    [0, 1],
    [5000, 0],
    ])

#plottimes: 11 evenly spaced points in time !only works with dt=0.5 or dt=1.0
plottimes = np.arange(0,TIME+dt,TIME/10)

#FOR LOOP, timetracker tt
tt = 0

for j in range(0, cols):
    # add load at time t
    for l, (time, load) in enumerate(tl):
        if tt == time:
            A[zero] += load

    #plot times of interest
    #plot von A, direkt GGÜ dzsum (dz ändert sich so immer genau richtig)
    if any(tt == plottimes):
        plt.plot(A[:],-dzsum)

    #iteration zeitvektoren:CALCULATING NEXT TIME STEP
    #Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
    B[zero] = fv[zero] * (f1[zero]*A[up] - 2*A[zero] + f2[zero]*A[lo]) + A[zero]
    A = B.copy()
    #timetracker: tt hat immer die Einheit der aktuellen Zeit in der Iteration (-> brauchbar für Zeiten des plots, und variable Lasten)
    tt += dt

#kontrolle der faktoren !<0.5
print('factors for each layer !<0.5')
ass.prnt_factors()

plt.show()