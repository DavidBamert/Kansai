"""
in: layer vectors 0-n
out: stacked/assembled vectors, ready for model iteration
"""

import numpy as np


class Assembly:
    def __init__(self, layerlist, dt, drainage, ob, yw):
        self.layerlist = layerlist
        self.dt = dt
        self.drainage = drainage
        self.ob = ob
        self.yw = yw

    def get_drainvect(self):
        drainvect = []
        for i in self.drainage:
            Lnew = self.layerlist[0:i]
            ssnew = Assembly(Lnew, self.dt, [], self.ob, self.yw)
            rownew = len(ssnew.get_dz())-1
            drainvect.append(rownew)
        return np.array(drainvect)

    def get_dzsum(self):
        array = np.empty((1,), float)
        hlow = 0
        for layer in self.layerlist:
            b = layer.get_dzsummed() + hlow
            array = np.delete(array, -1, 0)
            array = np.concatenate((array, b))
            hlow = layer.hlow
        return array

    def get_dz(self):
        array = np.empty((1,), float)
        for layer in self.layerlist:
            b = layer.get_dzvect()
            array = np.delete(array, -1, 0)
            array = np.concatenate((array, b))
        return array

    def get_k(self):
        array = np.empty((1,), float)
        for layer in self.layerlist:
            b = layer.get_kvect()
            array = np.delete(array, -1, 0)
            array = np.concatenate((array, b))
        return array

    def get_Me0_lin(self):
        array = np.empty((1,), float)
        for layer in self.layerlist:
            b = layer.get_Me0vect()
            array = np.delete(array, -1, 0)
            array = np.concatenate((array, b))
        return array

    def get_cv0_lin(self):
        array = np.empty((1,), float)
        for layer in self.layerlist:
            b = layer.get_cv0vect_lin()  # get constant cv for each different layer (linear calculation)
            array = np.delete(array, -1, 0)
            array = np.concatenate((array, b))
        return array

    # variable expansion for non-linear calculation
    def get_Cc(self):
        array = np.empty((1,), float)
        for layer in self.layerlist:
            b = layer.get_Cc()
            array = np.delete(array, -1, 0)
            array = np.concatenate((array, b))
        return array

    def get_e0(self):
        array = np.empty((1,), float)
        for layer in self.layerlist:
            b = layer.get_e0()
            array = np.delete(array, -1, 0)
            array = np.concatenate((array, b))
        return array

    def get_effsigma0(self):
        array = np.empty((1,), float)
        effsigma0low = 0
        for layer in self.layerlist:
            b = layer.get_effsigma0() + effsigma0low
            array = np.delete(array, -1, 0)
            array = np.concatenate((array, b))
            effsigma0low = array[-1]
        array += self.ob
        return array

    def get_Me0_nonlin(self):
        Me0 = np.log(10) * (1 + self.get_e0()) * self.get_effsigma0() / self.get_Cc()
        Me0[0] = Me0[1] / 2  # adjust the most upper Me, because it must not be 0
        return Me0

    # get factors for FDM
    def get_slices(self):
        # up = i-1; zero = i; lo = i+1
        rows = len(self.get_dz())
        up = slice(0, rows - 2)
        zero = slice(1, rows - 1)
        lo = slice(2, rows)
        return up, zero, lo

    def get_factors0(self, non_linear):
        # needed vectors
        dz = self.get_dz()
        k = self.get_k()
        if non_linear:
            cv0 = self.get_k() * self.get_Me0_nonlin() / self.yw  # non-linear calculation
        if not non_linear:
            cv0 = self.get_cv0_lin()  # linear calculation
        # length of vectors
        rows = len(self.get_dz())
        # up = i-1; zero = i; lo = i+1
        up, zero, lo = self.get_slices()
        # initialize
        fv = np.zeros((rows,))
        f1 = np.zeros((rows,))
        f2 = np.zeros((rows,))
        # factor vectors according to formula from LHAP page 66
        fv[zero] = np.array(((1 + k[up] / k[zero]) / (1 + cv0[zero] * k[up] / cv0[up] / k[zero])) * cv0[zero] * self.dt / dz[up] ** 2)
        f1[zero] = np.array(2 * k[up] / (k[zero] + k[up]))
        f2[zero] = np.array(2 * k[zero] / (k[zero] + k[up]))
        # complete first and last entry
        fv[0], fv[-1] = fv[1], fv[-2]
        f1[0], f1[-1] = f1[1], f1[-2]
        f2[0], f2[-1] = f2[1], f2[-2]
        return fv, f1, f2

    def get_mfact0_lin(self):
        mfact0 = []
        for layer in self.layerlist:
            factor = layer.cv0() * self.dt / layer.dzcorr() ** 2
            mfact0.append(factor)
        return np.array(mfact0)

    def get_mfact0_nonlin(self):
        mfact0 = []
        cv0 = self.get_k() * self.get_Me0_nonlin() / self.yw
        factors = cv0 * self.dt / self.get_dz() ** 2
        maxfact = max(factors)
        mfact0.append(maxfact)
        for i, fact in enumerate(factors):
            if fact == maxfact:
                maxdepth = self.get_dzsum()[i]

        return np.array(mfact0), maxdepth
