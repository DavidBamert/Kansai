""""
INPUT
damit dz kleiner gewählt werden kann, muss dt kleiner gewählt werden
"""

import timee as tm
import assembly as am
import layer as lm
import model as mm
import numpy as np


#discretizazion (dont use dt=0.3, for numerical noise reasons)
dz = 0.05
dt = 0.05

#Timeperiod
T = 100

#layers
L = [lm.Layer(0, 2, 0.3, 0.3, dz)
     ]

#add 1 dz to the last layers vector
L[-1].hlow += L[-1].dz

#boundry conditions [upper, lower] 0 drained, 1 undrained
bcs = [0, 0]

# loads in time tl = np.array([[time,load], ... ]) Matrix kann beliebig erweitert werden. Eintrag [0,1] kann IC ersetzen.
tl = np.array([
    [0, 1],
    [1000, 0]
])

#number of graphs
graphs = 11

#create assembly and timee object
ss = am.Assembly(L, dt)
tt = tm.Time(T, dt)

#solve the model using FDM
model = mm.Model(bcs, tl, ss, tt, graphs)
model.get_plot()


#assert all(tl[0] % dt == 0), 'tl[0] funktioniert nur, wenn von dt ohne rest geteilt'


"""
TODO:

Aufräumen

Plot verbessern:
    interpolieren
    A.reshape macht iteration evtl. langsamer
    
Drainage innerhalb der schichten: 
    model line 57 lösungsbeispiel
    A[h] = 0 
"""