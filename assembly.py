
"""
in: layervectors 0-n
out: stacked/assembled vectors, ready for model iteration
"""
import numpy as np
import layer as lm


class Assembly:
    def __init__(self, layerlist, dt, drainage):
        self.layerlist = layerlist
        self.dt = dt
        self.drainage = drainage

    def get_drainvect(self):
        drainvect = []
        for i in self.drainage:
            Lnew = self.layerlist[0:i]
            ssnew = Assembly(Lnew, self.dt, [])
            rownew = len(ssnew.get_dz())
            drainvect.append(rownew)

        return np.array(drainvect)



# dzsum hilft, die kurven unverzerrt zu plotten
    def get_dzsum(self):
        array = np.empty((0,), float)
        hlow = 0
        for layer in self.layerlist:
            b = layer.get_dzsummed() + hlow
            array = np.concatenate((array, b))
            hlow = layer.hlow
        return array

# vektoren, für iteration
    def get_dz(self):
        array = np.empty((0,), float)
        for layer in self.layerlist:
            b = layer.get_dzvect()
            array = np.concatenate((array, b))
        return array

    def get_k(self):
        array = np.empty((0,), float)
        for layer in self.layerlist:
            b = layer.get_kvect()
            array = np.concatenate((array, b))
        return array

    def get_cv(self):
        array = np.empty((0,), float)
        for layer in self.layerlist:
            b = layer.get_cvvect()
            array = np.concatenate((array, b))
        return array

# faktor vektoren und slices, für iteration

    def get_slices(self):
        # up = i-1; zero = i; lo = i+1
        rows = len(self.get_dz())

        up = slice(0, rows - 2)
        zero = slice(1, rows - 1)
        lo = slice(2, rows)
        return up, zero, lo

    def get_factors(self):
        #needed vectors
        dz = self.get_dz()
        k = self.get_k()
        cv = self.get_cv()

        #length of vectors
        rows = len(self.get_dz())
        #up = i-1; zero = i; lo = i+1
        up, zero, lo = self.get_slices()
        #initialize
        fv = np.zeros((rows,))
        f1 = np.zeros((rows,))
        f2 = np.zeros((rows,))
        #factor vectors according to formula s.66
        fv[zero] = np.array(((1 + k[up] / k[zero]) / (1 + cv[zero] * k[up] / cv[up] / k[zero])) * cv[zero] * self.dt / dz[up] ** 2)
        f1[zero] = np.array(2 * k[up] / (k[zero] + k[up]))
        f2[zero] = np.array(2 * k[zero] / (k[zero] + k[up]))
        #complete first and last entry
        fv[0], fv[-1] = fv[1], fv[-2]
        f1[0], f1[-1] = f1[1], f1[-2]
        f2[0], f2[-1] = f2[1], f2[-2]

        return fv, f1, f2

    # methode zum print aller factors der layers (!<0.5)
    def get_mfact(self):
        mfact = []
        for layer in self.layerlist:
            factor = layer.cv() * self.dt / layer.dzcorr() ** 2
            mfact.append(factor)
        return np.array(mfact)

# methode zum print aller vektoren, für kontrolle
    def prnt_vect(self):
        dzsum = self.get_dzsum()
        dz = self.get_dz()
        k = self.get_k()
        cv = self.get_cv()
        print(dzsum)
        print(dz)
        print(k)
        print(cv)
