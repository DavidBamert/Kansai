
import timee as tm
import assembly as am
import layer as lm
import numpy as np
import matplotlib.pyplot as plt

"""
wie kann ich am.Assembly nur 1 mal aufrufen, und nicht in jeder methode nochmal?
"""


class Model:
    def __init__(self, bcs, tl, ss, tt, amount):
        self.bcs = bcs  #1mal
        self.tl = tl    #ja
        self.ss = ss    #ja
        self.tt = tt    #ja
        self.amount = amount

    def get_fixeddata(self):

        #create A
        rows = len(self.ss.get_dz())
        A = np.zeros((rows,))

        #get factors and slices <-bcs
        factors_slices = self.ss.get_factors_slices(self.bcs)
        fv, f1, f2, up, zero, lo = factors_slices[:]
        #get time discretization and plot
        cols = self.tt.get_cols()
        plottimes, plotmatrix, timelegend = self.tt.get_plotmatrix(rows, self.amount)

        #ttrack = timetracker
        ttrack = 0
        #i = column number of plotmatrix
        i = 0

        return rows, A, fv, f1, f2, up, zero, lo, cols, plottimes, plotmatrix, timelegend, ttrack, i

    def get_variabledata(self):
        blubb = 1
        return blubb

    def solve(self):
            
        rows, A, fv, f1, f2, up, zero, lo, cols, plottimes, plotmatrix, timelegend, ttrack, i = self.get_fixeddata()

        # FOR LOOP,
        for j in range(0, cols):
            # add load at time t
            for l, (time, load) in enumerate(self.tl):
                if ttrack == time:
                    A[zero] += load

            # relevante Zeiten in matrix speichern
            #damit alle dt funktionieren: add 'if i<len(timelegend)''
            if ttrack >= plottimes[i]:
                plotmatrix[:, [i]] = np.reshape(A,(rows,1))
                timelegend[[i]] = ttrack
                i += 1

            # iteration zeitvektoren:CALCULATING NEXT TIME STEP
            # Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
            A[zero] = fv[zero] * (f1[zero] * A[up] - 2 * A[zero] + f2[zero] * A[lo]) + A[zero]
        
            # timetracker: tt hat immer die Einheit der aktuellen Zeit in der Iteration (-> brauchbar für Zeiten des plots, und variable Lasten)
            ttrack += self.tt.dt
            ttrack = round(ttrack, 3)   #dies macht dt robuster (eliminates numerical noise which causes problems)

        return plotmatrix, timelegend

    def get_plot(self):

        plotmatrix, timelegend = self.solve()

        # kontrolle der faktoren !<0.5
        print('factors for each layer !<0.5')
        self.ss.prnt_factors()

        # plot erstellen
        #dzsum helps plotting without distortion
        dzsum = self.ss.get_dzsum()
        plt.plot(plotmatrix[:], -dzsum, label = timelegend)
        plt.legend()
        plt.show()