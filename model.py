
import timee as tm
import assembly as am
import layer as lm
import numpy as np
import matplotlib.pyplot as plt

class Model:
    def __init__(self, dz, dt, T, L, bcs, tl):
        self.dz = dz
        self.dt = dt
        self.T = T
        self.L = L
        self.bcs = bcs
        self.tl = tl

    def solve(self):
            
        # use Assembly method and timee method
        ss = am.Assembly(self.L, self.dt)
        tt = tm.Time(self.T, self.dt)
        #get factors and slices <-bcs
        factors_slices = ss.get_factors_slices(self.bcs)
        fv, f1, f2, up, zero, lo = factors_slices[:]
        # time discretization
        cols = tt.get_cols()

        # IC, create A
        rows = len(ss.get_dz())
        A = np.zeros((rows,))

        # plottimes: 11 evenly spaced points in time !only works with dt=0.5 or dt=1.0
        plottimes = np.arange(0, self.T + self.dt, self.T/10)
        plotmatrix = np.zeros((rows,len(plottimes)))
        lege = np.zeros((len(plottimes),1))


        # FOR LOOP, timetracker tt, i for the plotmatrix columnnumber
        tt = 0
        i = 0

        for j in range(0, cols):
            # add load at time t
            for l, (time, load) in enumerate(self.tl):
                if tt == time:
                    A[zero] += load
        
            # relevante Zeiten in matrix speichern
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

        # plot erstellen
        #dzsum helps plotting without distortion
        dzsum = ss.get_dzsum()
        plt.plot(plotmatrix[:], -dzsum, label=lege)
        plt.legend()

        plt.show()
