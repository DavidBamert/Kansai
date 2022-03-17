# -*- coding: utf-8 -*-
"""
in: layer 0-n
out: assembled vectors, ready for iteration 
"""
import numpy as np
import layermethod as lm

#für beliebige Erweiterung: in: vektorliste; out: zusammengesetzte vektoren. 
#0:upper, 1:lower boundary (gemessen von oben)

class Assembly:
    def __init__(self, layerlist, dt):
        self.layerlist = layerlist
        self.dt = dt
        #dzsum hilft, die kurven unverzerrt zu plotten
    def get_dzsum(self):
        array = np.empty((0,),float)
        hlow = 0
        for L, (layer) in enumerate(self.layerlist):
            b = layer.get_dzsummed() + hlow
            array = np.concatenate((array , b))
            hlow = layer.hlow
        return array
        #folgend die vektoren, die für die iteration gebraucht werden
    def get_dz(self):
        array = np.empty((0,),float)
        for L, (layer) in enumerate(self.layerlist):
            b = layer.get_dzvect()
            array = np.concatenate((array , b))
        return array

    def get_k(self):
        array = np.empty((0,),float)
        for L, (layer) in enumerate(self.layerlist):
            b = layer.get_kvect()
            array = np.concatenate((array , b))
        return array

    def get_cv(self):
        array = np.empty((0,),float)
        for L, (layer) in enumerate(self.layerlist):
            b = layer.get_cvvect()
            array = np.concatenate((array , b))
        return array
        #methode zum print aller vektoren, für kontrolle
    def prnt_vect(self):
        dzsum= self.get_dzsum()
        dz = self.get_dz()
        k = self.get_k()
        cv = self.get_cv()
        print(dzsum)
        print(dz)
        print(k)
        print(cv)
        #methode zum print aller factors der layers
    def prnt_factors(self):
        for L, (layer) in enumerate(self.layerlist):
            factor= layer.cv()*self.dt/layer.dzcorr()**2
            print(round(factor,4))

