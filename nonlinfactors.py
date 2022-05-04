import timee as tm
import assembly as am
import layer as lm
import numpy as np


class eff_stress:
    def __init__(self, sigmaeff0):
        self.sigmaeff = sigmaeff0  #1mal

    def add_stress(self, du):
        self.sigmaeff += du

    def get_stress_eff(self):
        return np.array(self.sigmaeff)

class eff_e:
    def __init__(self, e0):
        self.e = e0

    def add_de(self, de):
        self.e += de

    def get_e_eff(self):
        return np.array(self.e)

"""
class dz:
    def __init__(self, dz):
        self.dz = dz

    def calculate strains(self, Me, delsigma):

    def calculate_new_dzvect(self):

    def get_dzvect(self):
        return np.array(dz)

"""