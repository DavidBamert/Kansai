"""
in: Layer
out: vectors with layerproperties
"""
# hup und hlow sind gemessen von oben (positiv)
# k: Durchlässigkeit kgross->sehr durchlÃ¤ssig
import math

import numpy as np
import input as ipt
import model as mdl

sigma_eff = ipt.tl #in order to make Me pressure dependend

yw = 10


class Layer:
    def __init__(self, hup, hlow, k, me_init, dz, c_c, e_0, me_case="constant"): #PD = pressureDependent
        self.hup = hup
        self.hlow = hlow
        self.me_case = me_case
        self.k = k
        self.me_init = me_init
        self.dz = dz
        self.c_c = c_c
        self.e_0 = e_0

#faktoren
    def height(self):
        height = self.hlow - self.hup
        return height

    def cv(self):
        cv = self.k * self.get_me() / yw
        return cv

    def get_me(self):
        i = mdl.get_fixeddata().row
        j = mdl.get_fixeddata().cols
        m_e = [i][j]

        if self.me_case == "constant": #me is constant
            self.m_e[i-1][j] = self.me_init
        else:   #me is pressure dependent
            if sigma_eff[i][j] <= sigma_eff[i][0] + 1 or sigma_eff[i][j] <= 0 or sigma_eff[i][0] <= 0: #to avoid a numerical error since potentially dividing by zero
                self.m_e[i-1][j] = (1+self.e_0)/self.c_c*(sigma_eff[i-1][0]+sigma_eff[i][0])/2
            else:
                self.m_e[i-1][j] = (1+self.e_0)/self.c_c*((sigma_eff[i-1][j]+sigma_eff[i][j])/2-(sigma_eff[i-1][0]+sigma_eff[i][0])/2)/(math.log10((sigma_eff[i-1][j]+sigma_eff[i][j])/2)-math.log10((sigma_eff[i-1][0]+sigma_eff[i][j])/2))
        return self.m_e

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
