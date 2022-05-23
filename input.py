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
bot = True
# non-linearity
nonlin = False
# second order strains
scnd = False

# discretizazion (dont use dt=0.3, for numerical noise reasons)
dz = 1
dt = 500

# time period
Tyears = 100
Tday = 365 * Tyears
T = 86400 * Tday

# layers (hup, hlow, k, Me, dz, gamma, Cc, e0, yw)
L = [lm.Layer(0, 12, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(12, 16, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(16, 28, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(28, 32, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(32, 44, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(44, 48, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(48, 60, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(60, 64, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(64, 76, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(76, 80, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(80, 92, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(92, 96, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(96, 108, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(108, 112, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(112, 124, 1e-9, 1670, dz, 8, 0.6, 1.5, yw),
     lm.Layer(124, 128, 1e-7, 1670, dz, 8, 0.1, 1.5, yw),
     lm.Layer(128, 140, 1e-9, 1670, dz, 8, 0.6, 1.5, yw)
     ]

# drainage inside the Layerassembly [1, 2, 3,....] (not more than layers-1 and >0)
drainage = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
dp = 0  # could be the waterpressure of a injection 'drainagepressure'
assert all(np.array(drainage) < len(L)) and all(np.array(drainage) > 0), 'more drainages than Layers-1'

# overburden pressure upon the modeled layers
ob = 200

# loads in time tl = np.array([[time,load], ... ])
# the matrix can be freely expanded, the entry [0,1] can replace the initial conditions
tl = np.array([
    [0, 403]
    ])

# number of graphs
graphs = 101  # number of exact solution u-vectors for U, Settlement s and the interpolation function

# create assembly and timee object
ss = am.Assembly(L, dt, drainage, ob, yw)
tt = tm.Time(T, dt)

# solve the model using FDM
model = mm.Model(tl, ss, tt, graphs, dp, yw)
solution = model.solve(top_drained=top, bot_drained=bot, non_linear=nonlin, sec_order_strains=scnd)
press = solution.plot_pressures(np.linspace(0, T, 11))
#solution.plot_pressures(np.linspace(0, T, 10), np.linspace(10, 20, 50)) # example of Urias
solution.plot_U()  # reference value 'U=1' is U(t=0)
settle = solution.plot_settlement(tl)

print('end')
