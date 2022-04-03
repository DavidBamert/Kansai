
import timee as tm
import assembly as am
import layer as lm
import numpy as np
import matplotlib.pyplot as plt

"""
wie kann ich am.Assembly nur 1 mal aufrufen, und nicht in jeder methode nochmal?
"""


class Model:
    def __init__(self, dz, dt, T, L, bcs, tl):
        self.dz = dz
        self.dt = dt
        self.T = T
        self.L = L
        self.bcs = bcs
        self.tl = tl

    def get_data(self):
        # use Assembly method and timee method
        ss = am.Assembly(self.L, self.dt)
        tt = tm.Time(self.T, self.dt)

        #create A
        rows = len(ss.get_dz())
        A = np.zeros((rows,))

        #get factors and slices <-bcs
        factors_slices = ss.get_factors_slices(self.bcs)
        fv, f1, f2, up, zero, lo = factors_slices[:]
        # time discretization
        cols = tt.get_cols()
        plottimes, plotmatrix, timelegend = tt.get_plotmatrix(rows)

        #tt = timetracker
        tt = 0
        #i = column number of plotmatrix
        i = 0

        return rows, A, fv, f1, f2, up, zero, lo, cols, plottimes, plotmatrix, timelegend, tt, i

    def solve(self):
            
        rows, A, fv, f1, f2, up, zero, lo, cols, plottimes, plotmatrix, timelegend, tt, i = self.get_data()

        # FOR LOOP,
        for j in range(0, cols):
            # add load at time t
            for l, (time, load) in enumerate(self.tl):
                if tt == time:
                    A[zero] += load

            # relevante Zeiten in matrix speichern
            #damit alle dt funktionieren: add 'if i<len(timelegend)''
            if tt >= plottimes[i]:
                plotmatrix[:, [i]] = np.reshape(A,(rows,1))
                timelegend[[i]] = tt
                i += 1

            # iteration zeitvektoren:CALCULATING NEXT TIME STEP
            # Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
            A[zero] = fv[zero] * (f1[zero] * A[up] - 2 * A[zero] + f2[zero] * A[lo]) + A[zero]
        
            # timetracker: tt hat immer die Einheit der aktuellen Zeit in der Iteration (-> brauchbar für Zeiten des plots, und variable Lasten)
            tt += self.dt
            tt = round(tt, 3)   #dies macht dt robuster (eliminates numerical noise which causes problems)

        return plotmatrix, timelegend

    def get_plot(self):

        ss = am.Assembly(self.L, self.dt)

        plotmatrix, timelegend = self.solve()

        # kontrolle der faktoren !<0.5
        print('factors for each layer !<0.5')
        ss.prnt_factors()

        # plot erstellen
        #dzsum helps plotting without distortion
        dzsum = ss.get_dzsum()
        plt.plot(plotmatrix[:], -dzsum, label = timelegend)
        plt.legend()
        plt.show()