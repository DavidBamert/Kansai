"""
in: Layer (hup und hlow are measured from above (positive))
out: vectors with layer properties
"""

import numpy as np


class Layer:
    def __init__(self, hup, hlow, k, me, dz, gamma, Cc, e0, yw):
        self.hup = hup
        self.hlow = hlow
        self.k = k
        self.me = me
        self.dz = dz
        self.gamma = gamma
        self.Cc = Cc
        self.e0 = e0
        self.yw = yw

    def cv(self):
        cv = self.k * self.me / self.yw
        return cv

# vectors
    # build vector dz, with an increment that is approximately = dz (so that h can be correctly represented)
    # is necessary for the plot, because dz differs from layer to layer
    def dzcorr(self):
        rows = int((self.hlow - self.hup) / self.dz)
        dzcorr = (self.hlow - self.hup) / rows
        return dzcorr

    def get_dzsummed(self):
        dzsumvect = np.arange(0, self.hlow - self.hup + self.dzcorr(), self.dzcorr())
        return dzsumvect

    def get_dzvect(self):
        dzvect = self.get_dzsummed()
        dzvect[:] = self.dzcorr()
        return dzvect

    def get_kvect(self):
        kvect = self.get_dzsummed()
        kvect[:] = self.k
        return kvect

    def get_cvvect(self):
        cv = self.k * self.me / self.yw
        cvvect = self.get_dzsummed()
        cvvect[:] = cv
        return cvvect

# variable expansion
    def get_effsigma(self):
        effsigma = self.gamma * self.get_dzsummed()
        return effsigma

    def get_Cc(self):
        Ccvect = self.get_dzsummed()
        Ccvect[:] = self.Cc
        return Ccvect

    def get_e0(self):
        e0vect = self.get_dzsummed()
        e0vect[:] = self.e0
        return e0vect
