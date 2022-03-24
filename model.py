# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 13:01:07 2022
damit dz kleiner gewählt werden kann, muss dt kleiner gewählt werden
"""
import assembly as am
import layer as lm
import numpy as np
import matplotlib.pyplot as plt

class Model:
    def __init__(self, dz, dt, TIME, L, bcs, tl):
        self.dz = dz
        self.dt = dt
        self.TIME = TIME
        self.L = L
        self.bcs = bcs
        self.tl =tl

    def solve(self):
            
        # throw L into Assembly method
        ss = am.Assembly(self.L, self.dt)
        
        # get all the iteration vektors
        dz = ss.get_dz()
        k = ss.get_k()
        cv = ss.get_cv()
        dzsum = ss.get_dzsum()
        
        # time discretization
        t = np.arange(0, self.TIME + self.dt, self.dt)
        cols = len(t)
        
        # IC, create B and A
        rows = len(dz)
        A = np.zeros((rows,))
        
        # BCs
        A[0] = self.bcs[0]
        A[-1] = self.bcs[1]
        
        """
        ITERATION FORWARD-DIFFERENCE METHOD !factor must be <0.5; 
        up = i-1
        zero = i
        lo = i+1
        """
        
        up = slice(0, rows - 2)
        zero = slice(1, rows - 1)
        lo = slice(2, rows)
        
        # factor vectors must be the same length as the dzsum vector
        
        fv = np.zeros((rows,))
        f1 = np.zeros((rows,))
        f2 = np.zeros((rows,))
        
        # dzn[:] = np.array(dz)
        fv[zero] = np.array(((1 + k[zero] / k[lo]) / (1 + cv[lo] * k[zero] / cv[zero] / k[lo])) * cv[lo] * self.dt / dz[zero] ** 2)
        f1[zero] = np.array(2 * k[zero] / (k[lo] + k[zero]))
        f2[zero] = np.array(2 * k[lo] / (k[lo] + k[zero]))
        
        # complete the first and last entry
        fv[0], fv[-1] = fv[1], fv[-2]
        f1[0], f1[-1] = f1[1], f1[-2]
        f2[0], f2[-1] = f2[1], f2[-2]

        # plottimes: 11 evenly spaced points in time !only works with dt=0.5 or dt=1.0
        plottimes = np.arange(0, self.TIME + self.dt, self.TIME/10)
        plotmatrix = np.zeros((rows,len(plottimes)))
        lege = np.zeros((len(plottimes),1))
        
        # FOR LOOP, timetracker tt
        tt = 0
        i = 0
        
        for j in range(0, cols):
            # add load at time t
            for l, (time, load) in enumerate(self.tl):
                if tt == time:
                    A[zero] += load
        
            # plot times of interest
            # plot von A, direkt GGÜ dzsum (dz ändert sich so immer genau richtig)
            if any(plottimes == tt):
                plotmatrix[:, [i]] = np.reshape(A,(rows,1))
                lege[[i]] = tt
                i += 1
        
            # iteration zeitvektoren:CALCULATING NEXT TIME STEP
            # Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
            A[zero] = fv[zero] * (f1[zero] * A[up] - 2 * A[zero] + f2[zero] * A[lo]) + A[zero]
        
            # timetracker: tt hat immer die Einheit der aktuellen Zeit in der Iteration (-> brauchbar für Zeiten des plots, und variable Lasten)
            tt += self.dt
        
        # kontrolle der faktoren !<0.5
        print('factors for each layer !<0.5')
        ss.prnt_factors()
        
        plt.plot(plotmatrix[:], -dzsum, label=lege)
        plt.legend()
        
        plt.show()
