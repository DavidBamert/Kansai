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
dz = 0.1
dt = 2e4

#Timeperiod
Tday = 5e3
T = 86400 * Tday

#layers(self, hup, hlow, k, me, dz, gamma, Cc, e0, yw)
L = [lm.Layer(0, 12, 1e-9, 1670, dz, 10, 0.6, 1.5, yw)
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
uexact = 101           #number of exact solution u-vectors for U, Settlement s and the interpolation function

#create assembly and timee object
ss = am.Assembly(L, dt, drainage, yw)
tt = tm.Time(T, dt)

#solve the model using FDM
model = mm.Model(tl, ss, tt, uexact, dp, yw)
solution = model.solve(bot_drained=True)
press = solution.plot_pressures(np.linspace(0, T, 5))
#solution.plot_pressures(np.linspace(0, T, 10), np.linspace(10, 20, 50)) #beispiel von urias

solution.plot_U()            #Referenzwert 'U=1' ist U(t=0)
settle = solution.plot_settlement(tl)

print('end')

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