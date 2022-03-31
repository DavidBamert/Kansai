
import numpy as np
import layer as lm

class Time:
    def __init__(self, T, dt):
        self.T = T
        self.dt = dt

    def get_cols(self):
        t = np.arange(0, self.T + self.dt, self.dt)
        cols = len(t)
        return cols