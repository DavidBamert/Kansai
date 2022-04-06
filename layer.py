"""
in: Layer
out: vectors with layerproperties
"""
# hup und hlow sind gemessen von oben (positiv)
# k: Durchlässigkeit kgross->sehr durchlÃ¤ssig

import numpy as np

yw = 10


class Layer:
    def __init__(self, hup, hlow, k, me, dz):
        self.hup = hup
        self.hlow = hlow
        self.k = k
        self.me = me
        self.dz = dz

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

#zur kontrolle
    def prnt_vect(self):
        dzvect = self.get_dzvect()
        kvect = self.get_kvect()
        cvvect = self.get_cvvect()
        print(dzvect)
        print(kvect)
        print(cvvect)
