# -*- coding: utf-8 -*-
"""
INPUT
"""

import assembly as am
import layer as lm
import model as mm
import numpy as np


#discretizazion
dz = 0.5
dt = 0.5

#Timeperiod
TIME = 10000

#layers
L = [lm.Layer(0, 10, 0.2, 0.3, dz),
     lm.Layer(10, 20, 0.2, 0.3, dz),
     lm.Layer(20, 30, 0.2, 0.3, dz)
     ]

#boundry conditions [upper, lower] 0 drained, 1 undrained
bcs = [0, 0]

# loads in time tl = np.array([[time,load], ... ]) Matrix kann beliebig erweitert werden. Eintrag [0,1] kann IC ersetzen.
tl = np.array([
    [0, 1],
    [5000, 0],
])

model = mm.Model(dz,dt,TIME,L,bcs,tl)
model.solve()

"""
todo:
viel von modelmethod in assemblymetohd verschieben (Faktorvektoren, Plotcodes,...)
etwa 1000 zeitvektoren speichern, und für plot dann relevante raussuchen 

undrained cond: k[0] k[-1] = 0 setzen 



"""