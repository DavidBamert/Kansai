"""
in: -ganze Höhe, diskretisieren -> Z-vektor
    -übergangsbedingung = allgemeiner fall! nur diese Formel im code verwenden? 

    -einzelne Layers: hsup, hinf, <- in diesem intervall werden Cv und k definiert (in Vektor mit länge Z-vektor ->k-Vektor ->Cv Vektor)
        für Rechnung an stelle i fragt die Formel immer k1(Cv1) and stelle i und k2(Cv2) and stelle i+1 ab! 
            veränderbares dz ist auch so lösbar
    
deltaZ verschieden in layers: 
    +schichtgrenzen sind genau
    -dz ist variabel in formel 
"""

import numpy as np

yw = 10
dz = 0.5

#hup und hlow sind gemessen von oben (positiv)
#k: Durchlässigkeit kgross->sehr durchlässig
class layer:
    def __init__(self, hup, hlow, k, me):
        self.hup = hup
        self.hlow = hlow
        self.k = k
        self.me= me 
        
        
    def height(self):
        height = self.hlow-self.hup
        return height
    
    def cv(self):
        cv = self.k*self.me/yw
        return cv
    
        #vector dz erstellen, mit inkrement das ungefähr = dz ist. (so dass h genau abgebildet werden kann)       
    def dzcorr(self):
        rows = int((self.hlow-self.hup)/dz)
        dzcorr = (self.hlow-self.hup)/rows
        return dzcorr

        #ist nötig für den plot, deil dz von layer zu layer anders ist.
    def get_dzsummed(self):
        dzsumvect = np.arange(0,self.hlow-self.hup+self.dzcorr(),self.dzcorr())
        return dzsumvect
        """
        Vektoren für die Iteration
        """
    def get_dzvect(self):
        dzvect = self.get_dzsummed()
        dzvect[:] = self.dzcorr()
        return dzvect
    
    def get_kvect(self):
        kvect = self.get_dzsummed()
        kvect[:] = self.k
        return kvect
    
    def get_cvvect(self):
        cv = self.k*self.me/yw
        cvvect = self.get_dzsummed()
        cvvect[:] = cv
        return cvvect 
        
        """
        zur kontrolle
        """
    def prnt_vect(self):
        dzvect = self.get_dzvect()
        kvect = self.get_kvect()
        cvvect = self.get_cvvect()
        print(dzvect)
        print(kvect)
        print(cvvect)
    
#yolo = layer(1,10,1,1)
#yolo.prnt_vect()