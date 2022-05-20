"""
INPUT
"""

import timee as tm
import assembly as am
import layer as lm
import model as mm
import numpy as np

yw = 10
# drainage
top = True
bot = False
# second order strains
scnd = False

# discretizazion (dont use dt=0.3, for numerical noise reasons)
dz = 1
dt = 500

# time period
Tday = 18250  # =50 years
T = 86400 * Tday

# layers(self, hup, hlow, k, me, dz, gamma, Cc, e0, yw)
L = [lm.Layer(0, 21, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(21, 25, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(25, 37, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(37, 41, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(41, 53, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(53, 57, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(57, 69, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(69, 73, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(73, 85, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(85, 89, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(89, 101, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(101, 105, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(105, 117, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(117, 121, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(121, 133, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(133, 137, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(137, 149, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(149, 153, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(153, 165, 1e-9, 1670, dz, 8, 0.6, 1.5, yw)
     ]

# drainage inside the Layerassembly [1, 2, 3,....] (not more than layers-1 and >0)
drainage = []
dp = 0  # could be the waterpressure of a injection 'drainagepressure'
assert all(np.array(drainage) < len(L)) and all(np.array(drainage) > 0), 'more drainages than Layers-1'

# loads in time tl = np.array([[time,load], ... ])
# the matrix can be freely expanded, the entry [0,1] can replace the initial conditions
tl = np.array([
    [0, 403]
    ])

# number of graphs
uexact = 101  # number of exact solution u-vectors for U, Settlement s and the interpolation function

# create assembly and timee object
ss = am.Assembly(L, dt, drainage, yw)
tt = tm.Time(T, dt)

# solve the model using FDM
model = mm.Model(tl, ss, tt, uexact, dp, yw)
solution = model.solve(top_drained=top, bot_drained=bot, sec_order_strains=scnd)
press = solution.plot_pressures(np.linspace(0, T, 5))
#solution.plot_pressures(np.linspace(0, T, 10), np.linspace(10, 20, 50))  # example of Urias

solution.plot_U()  # reference value 'U=1' is U(t=0)
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