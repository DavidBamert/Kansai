"""
INPUT
To choose smaller dz, dt must be decreased much more
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
# non-linearity
nonlin = False
# second order strains
scnd = False

# discretizazion (dont use dt=0.3, for numerical noise reasons)
dz = 1
dt = 500

# time period
Tyears = 1
Tday = 365 * Tyears
T = 86400 * Tday

# layers (hup, hlow, k, me, dz, gamma, Cc, e0, yw)
L = [lm.Layer(0, 17.5, 1e-8, 1670, dz, 5, 0.6, 1.5, yw),
     lm.Layer(17.5, 30, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(30, 32.5, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(32.5, 43, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(43, 50.5, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(50.5, 55, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(55, 57, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(57, 61, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(61, 62.5, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(62.5, 82, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(82, 86.5, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(86.5, 94, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(94, 99, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(99, 103, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(103, 106, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(106, 111, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(111, 118, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(118, 120, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(120, 131, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(131, 157, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(157, 172, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     ]


# drainage inside the Layerassembly [1, 2, 3,....] (not more than layers-1 and >0)
drainage = []
dp = 0  # could be the waterpressure of a injection 'drainagepressure'
assert all(np.array(drainage) < len(L)) and all(np.array(drainage) > 0), 'more drainages than Layers-1'

# overburden pressure upon the modeled layers
ob = 0

# loads in time tl = np.array([[time,load], ... ])
# the matrix can be freely expanded, the entry [0,1] can replace the initial conditions
tl = np.array([
    [0, 40]
    ])

# number of graphs
uexact = 101  # number of exact solution u-vectors for U, Settlement s and the interpolation function

# create assembly and timee object
ss = am.Assembly(L, dt, drainage, ob, yw)
tt = tm.Time(T, dt)

# solve the model using FDM
model = mm.Model(tl, ss, tt, uexact, dp, yw)
solution = model.solve(top_drained=top, bot_drained=bot, non_linear=nonlin, sec_order_strains=scnd)
press = solution.plot_pressures(np.linspace(0, T, 100))
#solution.plot_pressures(np.linspace(0, T, 10), np.linspace(10, 20, 50))  # example of Urias
solution.plot_U()  # reference value 'U=1' is U(t=0)
settle = solution.plot_settlement(tl)

print('end')
