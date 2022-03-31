# -*- coding: utf-8 -*-
"""
INPUT
damit dz kleiner gewählt werden kann, muss dt kleiner gewählt werden
"""

import assembly as am
import layer as lm
import model as mm
import numpy as np


#discretizazion
dz = 0.5
dt = 0.5

#Timeperiod
T = 100000

#layers
L = [lm.Layer(0, 20, 0.2, 0.3, dz),
     lm.Layer(20, 40, 0.2, 0.3, dz),
     lm.Layer(40, 60, 0.2, 0.3, dz)
     ]

#add 1 dz to the last layers vector
L[-1].hlow += L[-1].dz

#boundry conditions [upper, lower] 0 drained, 1 undrained
bcs = [0, 0]

# loads in time tl = np.array([[time,load], ... ]) Matrix kann beliebig erweitert werden. Eintrag [0,1] kann IC ersetzen.
tl = np.array([
    [0, 1],
    [25000, 0],
])

model = mm.Model(dz, dt, T, L, bcs, tl)
model.solve()


"""
TODO:

Aufräumen:
    viel von modelmethod in assemblymetohd verschieben (Faktorvektoren, Plotcodes,...)

Plot verbessern:
    etwa n*1000 zeitvektoren speichern, und für plot die relevanten raussuchen. interpolieren (?)

undrained cond: 
    k[0] k[-1] = 0 setzen 
    verlust von 1dz, wenn wir es so machen. layer wird dünner
    bessere lösung nötig


"""