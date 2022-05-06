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
dz = 0.5
dt = 100000

#Timeperiod
T = 1e9

#layers(self, hup, hlow, k, me, dz, gamma, Cc, e0)
L = [lm.Layer(0, 10, 1e-9, 1700, dz, 10, 0.4, 0.9),
     lm.Layer(10, 20, 1e-9, 1700, dz, 20, 0.4, 0.9),
     lm.Layer(20, 30, 1e-9, 1700, dz, 12, 0.4, 0.9)
     ]

#drainage inside the Layerassembly [1, 2, 3,....] (not more than layers-1 and >0)
drainage = []
dp = 0 #could be the waterpressure of a injection 'drainagepressure'
assert all(np.array(drainage) < len(L)) and all(np.array(drainage) > 0), 'more drainages than Layers-1'

#boundry conditions [upper, lower] 0 drained, 1 undrained
bcs = [0, 0]
assert bcs == [0, 0] or bcs == [0, 1] or bcs == [1, 0] or bcs == [1, 1], 'check bcs'

# loads in time tl = np.array([[time,load], ... ]) Matrix kann beliebig erweitert werden. Eintrag [0,1] kann IC ersetzen.
tl = np.array([
    [0, 100],
    [0.5e9, 20]
    ])

#number of graphs
graphs = 1001           #number of exact solution vectors for U and the interpolation function

#create assembly and timee object
ss = am.Assembly(L, dt, drainage)
tt = tm.Time(T, dt)

#solve the model using FDM
model = mm.Model(bcs, tl, ss, tt, graphs, dp)
solution = model.solve(top_drained=True)
solution.plot_pressures(np.linspace(0, T, 11))
#solution.plot_pressures(np.linspace(0, T, 10), np.linspace(10, 20, 50)) #beispiel von urias

solution.plot_U()            #Referenzwert 'U=1' ist U(t=0)
solution.plot_settlement()

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
"""

print(ss.get_effsigma())
print(ss.get_e0())
print(ss.get_Cc())

print(ss.get_Me())
"""