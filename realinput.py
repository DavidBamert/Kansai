"""
INPUT
"""

import timee as tm
import assembly as am
import layer as lm
import model as mm
import numpy as np


#discretizazion (dont use dt=0.3, for numerical noise reasons)
dz = 0.5
dt = 1e4

#Timeperiod
T = 5e9

#layers
L = [lm.Layer(0, 20, 1e-9, 1670, dz, 11, 0.4, 1.0),
     lm.Layer(20, 40, 1e-9, 1670, dz, 11, 0.4, 1.0)
     ]

#add 1 dz to the last layers vector
L[-1].hlow += L[-1].dz

#drainage inside the Layerassembly [1, 2, 3,....] (not more than layers-1)
drainage = []
dp = 0 #could be the waterpressure of a injection 'drainagepressure'
assert all(np.array(drainage) < len(L)), 'more drainages than Layers-1'

#boundry conditions [upper, lower] 0 drained, 1 undrained
bcs = [0, 0]
assert bcs == [0, 0] or bcs == [0, 1] or bcs == [1, 0] or bcs == [1, 1], 'check bcs'

# loads in time tl = np.array([[time,load], ... ]) Matrix kann beliebig erweitert werden. Eintrag [0,1] kann IC ersetzen.
tl = np.array([
    [0, 1],
    [3e9, 0]
])

#number of graphs
graphs = 11

#create assembly and timee object
ss = am.Assembly(L, dt, drainage)
tt = tm.Time(T, dt)

model = mm.Model(bcs, tl, ss, tt, graphs, dp)
model.get_plot()

