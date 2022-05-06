
import timee as tm
import assembly as am
import layer as lm
import model as mm
import numpy as np


dz = 0.5
dt = 5e3  # seconds

# Time period
T = 1e9  # must be in days

# Layers (hup, hlow, k, me, dz, gamma, Cc, e0)
L = [lm.Layer(0, 12, 1e-9, 1700, dz, 12, 0.4, 0.9),
     lm.Layer(12, 16, 1e-7, 1700, dz, 9, 0.4, 0.9),
     lm.Layer(16, 28, 1e-9, 1700, dz, 12, 0.4, 0.9)
     ]

# Drainage inside the Layerassembly [1, 2, 3,....] (not more than layers-1 and >0)
drainage = []
dp = 0  # could be the waterpressure of a injection 'drainagepressure'
assert all(np.array(drainage) < len(L)) and all(np.array(drainage) > 0), 'more drainages than Layers-1'

# Boundary conditions [upper, lower] 0 drained, 1 undrained
bcs = [0, 0]
assert bcs == [0, 0] or bcs == [0, 1] or bcs == [1, 0] or bcs == [1, 1], 'check bcs'

# Loads in time tl = np.array([[time,load], ... ])
# Matrix kann beliebig erweitert werden. Eintrag [0,1] kann IC ersetzen.
tl = np.array([
    [0, 1],
    [2.5e9, 0]
])

# Number of graphs
graphs = 11  # number of exact solution vectors for U and the interpolation function

# Create assembly and timee object
ss = am.Assembly(L, dt, drainage)
tt = tm.Time(T, dt)  # also graphs?

ss.get_mfact()

# Solve the model using FDM
model = mm.Model(tl, ss, tt, graphs, dp)
solution = model.solve(bot_drained=True)
solution.plot_pressures(np.linspace(0, T, 11))
# solution.plot_pressures(np.linspace(0, T, 10), np.linspace(10, 20, 50))

solution.get_U()  # Referenzwert 'U=1' ist U(t=0)
solution.get_dzz()

print('end')
