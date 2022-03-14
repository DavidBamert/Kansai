# -*- coding: utf-8 -*-
"""
in: layer 0-n
out: assembled vectors, ready for iteration 
"""
import numpy as np
import layermethod as lm

#für beliebige Erweiterung: in: vektorliste; out: zusammengesetzte vektoren. 
#0:upper, 1:lower boundary (gemessen von oben)

layer0 = lm.layer(0,2,0.04,0.3)
layer1 = lm.layer(2,3,0.04,0.3)


class assembly:
    def __init__(self,vectors):
        self.vectors = vectors
        
        """
        bei dzsum ist wichtig, dass der Vektor keine Stelle wiederholt-> sichtbar als sprung im plot. Deshalb +layer1.height() +layer1.dzcorr()
        """
    def get_dzsum(self):
        a,b = layer0.get_dzsummed(),layer1.get_dzsummed()+layer0.height()+layer0.dzcorr()
        cv = np.concatenate((a,b))
        return cv        
    
    def get_dz(self):
        a,b = layer0.get_dzvect(),layer1.get_dzvect()
        dz = np.concatenate((a,b))
        return dz
            
    def get_k(self):
        a,b = layer0.get_kvect(),layer1.get_kvect()
        k = np.concatenate((a,b))
        return k
        
    def get_cv(self):
        a,b = layer0.get_cvvect(),layer1.get_cvvect()
        cv = np.concatenate((a,b))
        return cv
    
    def prnt_vect(self):
        dzsum= self.get_dzsum()
        dz = self.get_dz()
        k = self.get_k()
        cv = self.get_cv()
        print(dzsum)
        print(dz)
        print(k)
        print(cv)
        



"""
nur zur kontrolle
"""

dt=0.1
factor0= layer0.cv()*dt/layer0.dzcorr()**2
factor1= layer1.cv()*dt/layer1.dzcorr()**2
print('Faktoren für die verschiedenen Layers !<0.5')
print(factor0)
print(factor1)