"""
this class solves the time discretization and the plot matrix
"""


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

    #plottimes: save vector A as soon as tt>=plottime
    #plotmatrix: here the vector A is saved
    #timelegend: here the precise time is saved (this is not always equal to plottimes)
    def get_plotmatrix(self, rows, amount):

        plottimes = np.linspace(0, self.T, amount) #<- here the amount of saved vectors can be chosen
        plotmatrix = np.zeros((rows, len(plottimes)))
        timelegend = np.zeros((len(plottimes), 1))
        return plottimes, plotmatrix, timelegend