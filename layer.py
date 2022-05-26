"""
in: Layer (hup und hlow are measured from above (positive))
out: vectors with layer properties
"""

import numpy as np


class Layer:
    def __init__(self, hup, hlow, k, Me0, dz, gamma, Cc, e0, yw):
        self.hup = hup
        self.hlow = hlow
        self.k = k
        self.Me0 = Me0
        self.dz = dz
        self.gamma = gamma
        self.Cc = Cc
        self.e0 = e0
        self.yw = yw

    def cv0(self):
        cv0 = self.k * self.Me0 / self.yw  # cv calculated from user-defined Me
        return cv0

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

    def get_Me0vect(self):
        Me0vect = self.get_dzsummed()
        Me0vect[:] = self.Me0
        return Me0vect

    def get_cv0vect_lin(self):
        cv0 = self.k * self.Me0 / self.yw  # cv calculated from user-defined Me0
        cv0vect = self.get_dzsummed()
        cv0vect[:] = cv0
        return cv0vect  # for linear calculation

    # variable expansion for non-linear calculation
    def get_Cc(self):
        Ccvect = self.get_dzsummed()
        Ccvect[:] = self.Cc
        return Ccvect

    def get_e0(self):
        e0vect = self.get_dzsummed()
        e0vect[:] = self.e0
        return e0vect

    def get_effsigma0(self):
        effsigma0 = self.gamma * self.get_dzsummed()
        return effsigma0
