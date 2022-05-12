"""
INPUT
damit dz kleiner gewählt werden kann, muss dt kleiner gewählt werden
"""

import timee as tm
import assembly as am
import layer as lm
import model as mm
import numpy as np

yw = 10

#discretizazion (dont use dt=0.3, for numerical noise reasons)
dz = 0.5
dt = 1e4

#Timeperiod
T = 1.728e9

#layers(self, hup, hlow, k, me, dz, gamma, Cc, e0, yw)
L = [lm.Layer(0, 12, 1e-9, 1670, dz, 12, 0.6, 1.5, yw),
    lm.Layer(12, 16, 1e-7, 1670, 1, 10, 0.6, 1.5, yw),
    lm.Layer(16, 28, 1e-9, 1670, dz, 15, 0.6, 1.5, yw)
     ]

#drainage inside the Layerassembly [1, 2, 3,....] (not more than layers-1 and >0)
drainage = []
dp = 0 #could be the waterpressure of a injection 'drainagepressure'
assert all(np.array(drainage) < len(L)) and all(np.array(drainage) > 0), 'more drainages than Layers-1'

# loads in time tl = np.array([[time,load], ... ]) Matrix kann beliebig erweitert werden. Eintrag [0,1] kann IC ersetzen.
tl = np.array([
    [0, 100]
    ])

#number of graphs
graphs = 1001           #number of exact solution u-vectors for U, Settlement s and the interpolation function

#create assembly and timee object
ss = am.Assembly(L, dt, drainage, yw)
tt = tm.Time(T, dt)

#solve the model using FDM
model = mm.Model(tl, ss, tt, graphs, dp, yw)
solution = model.solve(bot_drained=True)
solution.plot_pressures(np.linspace(0, T, 11))
#solution.plot_pressures(np.linspace(0, T, 10), np.linspace(10, 20, 50)) #beispiel von urias

solution.plot_U()            #Referenzwert 'U=1' ist U(t=0)
solution.plot_settlement2(tl)



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