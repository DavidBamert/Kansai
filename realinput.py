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
nonlin = True
# second order strains
scnd = True

# discretizazion (dont use dt=0.3, for numerical noise reasons)
dz = 0.5
dt = 50

# time period
Tyears = 100
Tday = 1000 #365 * Tyears
T = 86400 * Tday

maxk = 1e-7

# layers (hup, hlow, k, me, dz, gamma, Cc, e0, yw)
L = [lm.Layer(0,    17.5,  8.44e-9, 1670, dz, 5.1, 0.69, 2.30, yw),
     lm.Layer(17.5, 30,       maxk, 1670, dz, 9.6, 0.03, 0.40, yw),
     lm.Layer(30,   32.5,  1.35e-9, 1670, dz, 6.7, 0.43, 1.19, yw),
     lm.Layer(32.5, 43,    8.7e-10, 1670, dz, 5.1, 0.72, 1.78, yw),
     lm.Layer(43,   50.5,     maxk, 1670, dz, 9.6, 0.03, 0.40, yw),
     lm.Layer(50.5, 55,   4.45e-10, 1670, dz, 6.5, 0.54, 1.31, yw),
     lm.Layer(55,   57,       maxk, 1670, dz, 9.6, 0.03, 0.40, yw),
     lm.Layer(57,   61,    8.2e-10, 1670, dz, 6.7, 0.54, 1.33, yw),
     lm.Layer(61,   62.5,     maxk, 1670, dz, 9.6, 0.03, 0.40, yw),
     lm.Layer(62.5, 82,    6.2e-10, 1670, dz, 6.5, 0.70, 1.42, yw),
     lm.Layer(82,   86.5,     maxk, 1670, dz, 9.6, 0.03, 0.40, yw),
     lm.Layer(86.5, 94,   4.55e-10, 1670, dz, 6.5, 0.72, 1.43, yw),
     lm.Layer(94,   99,       maxk, 1670, dz, 9.6, 0.03, 0.40, yw),
     lm.Layer(99,   103,   3.4e-10, 1670, dz, 5.8, 0.76, 1.51, yw),
     lm.Layer(103,  106,      maxk, 1670, dz, 9.6, 0.03, 0.40, yw),
     lm.Layer(106,  111,     6e-10, 1670, dz, 7.1, 0.61, 1.29, yw),
     lm.Layer(111,  118,      9e-8, 1670, dz, 9.6, 0.03, 0.40, yw),
     lm.Layer(118,  120,  3.43e-10, 1670, dz, 6.7, 0.69, 1.30, yw),
     lm.Layer(120,  131,     4e-10, 1670, dz, 9.6, 0.03, 0.40, yw),  # 19
     lm.Layer(131,  157,  3.69e-10, 1670, dz, 6.4, 0.69, 1.28, yw),
     lm.Layer(157,  172,  2.78e-10, 1670, dz, 6.4, 0.55, 0.84, yw)
     ]

# drainage inside the Layerassembly [1, 2, 3,....] (not more than layers-1 and >0)
drainage = [18, 19]
dp = 0  # could be the waterpressure of a injection 'drainagepressure'
assert all(np.array(drainage) < len(L)) and all(np.array(drainage) > 0), 'more drainages than Layers-1'

# overburden pressure upon the modeled layers
ob = 0

# loads in time tl = np.array([[time,load], ... ])
# the matrix can be freely expanded, the entry [0,1] can replace the initial conditions
ts = 86400  # seconds in one day
tl = np.array([
     [0,        30],
     [175*ts,   30],
     [655*ts,  115],
     [1075*ts,  40],
     [1195*ts, 140],
     [1255*ts,  60],
     [1855*ts,  45]
    ])

# number of graphs
uexact = 101  # number of exact solution u-vectors for U, Settlement s and the interpolation function

# create assembly and timee object
ss = am.Assembly(L, dt, drainage, ob, yw)
tt = tm.Time(T, dt)

# solve the model using FDM
model = mm.Model(tl, ss, tt, uexact, dp, yw)
solution = model.solve(top_drained=top, bot_drained=bot, non_linear=nonlin, sec_order_strains=scnd)
umatrix = solution.plot_pressures(np.linspace(0, T, 11))
#solution.plot_pressures(np.linspace(0, T, 10), np.linspace(10, 20, 50))  # example of Urias
solution.plot_U()  # reference value 'U=1' is U(t=0)
settlementvect = solution.plot_settlement(tl)

print('end')
