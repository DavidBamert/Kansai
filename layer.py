"""
in: Layer
out: vectors with layerproperties
"""
# hup und hlow sind gemessen von oben (positiv)
# k: Durchlässigkeit kgross->sehr durchlÃ¤ssig

import numpy as np

yw = 10


class Layer:
    def __init__(self, hup, hlow, k, me, dz, gamma, Cc, e0):
        self.hup = hup
        self.hlow = hlow
        self.k = k
        self.me = me
        self.dz = dz
        self.gamma = gamma
        self.Cc = Cc
        self.e0 = e0

#faktoren
    def height(self):
        height = self.hlow - self.hup
        return height

    def cv(self):
        cv = self.k * self.me / yw
        return cv

        # vector dz erstellen, mit inkrement das ungefÃ¤hr = dz ist. (so dass h genau abgebildet werden kann)
        # ist nÃ¶tig fÃ¼r den plot, deil dz von layer zu layer anders ist.

#vektoren
    def dzcorr(self):
        rows = int((self.hlow - self.hup) / self.dz)
        dzcorr = (self.hlow - self.hup) / rows
        return dzcorr

        #Vektoren für die Iteration

    def get_dzsummed(self):
        dzsumvect = np.arange(0, self.hlow - self.hup, self.dzcorr())
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
        cv = self.k * self.me / yw
        cvvect = self.get_dzsummed()
        cvvect[:] = cv
        return cvvect

#variable expansion
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

#zur kontrolle
    def prnt_vect(self):
        dzvect = self.get_dzvect()
        kvect = self.get_kvect()
        cvvect = self.get_cvvect()
        print(dzvect)
        print(kvect)
        print(cvvect)
