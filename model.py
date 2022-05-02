
import timee as tm
import assembly as am
import layer as lm
import numpy as np
import matplotlib.pyplot as plt
yw = 10

class Model:
    def __init__(self, bcs, tl, ss, tt, graphs, dp):
        self.bcs = bcs  #1mal
        self.tl = tl    #ja
        self.ss = ss    #ja
        self.tt = tt    #ja
        self.graphs = graphs
        self.dp = dp
    def get_fixeddata(self):    #this prepares the fix data for the iteration

        #create A
        rows = len(self.ss.get_dz())
        A = np.zeros((rows,))

        #get factors and slices <-bcs
        up, zero, lo = self.ss.get_slices()

        #get time discretization and plot
        cols = self.tt.get_cols()
        plottimes, plotmatrix, timelegend = self.tt.get_plotmatrix(rows, self.graphs)

        #ttrack = timetracker
        ttrack = 0
        #i = column number of plotmatrix
        i = 0

        return rows, A, up, zero, lo, cols, plottimes, plotmatrix, timelegend, ttrack, i

    def get_factorvectors(self): #this prepares (and in the future alters) the variable data for the iteration
        #change fi and dz
        fv, f1, f2 = self.ss.get_factors()
        return fv, f1, f2

    def get_updated_factorvectors(self, deltau, udisstot):


        #needed vectors
        Me = self.ss.get_Me()
        dz = self.ss.get_dz()
        k = self.ss.get_k()
        e0 = self.ss.get_e0()
        Cc = self.ss.get_Cc()
        #calculate sigma1 and sigma2
        sigma2 = self.ss.get_effsigma() + udisstot
        sigma1 = sigma2 - deltau
        #calculate new e
        j=0
        for sigma11,sigma22 in zip(sigma1,sigma2): #wie kann man hier fallunterscheiden? !!delta e0 aus jeder iteration müsste summiert sein!!
            if sigma22-sigma11 > 1e-7:
                e0[j] = e0[j] - Cc[j] * np.log10(sigma22/sigma11)
                j +=1
        #calculate new Me
        i=0
        for sigma11,sigma22 in zip(sigma1,sigma2): #problem with log <- fallunterscheidung für die Me berechnung
            if sigma22-sigma11 < 1e-7:
                Me[i] = np.log(10) * (1+e0[i]) / Cc[i] * sigma22
                i +=1
            else:
                Me[i] = (1 + e0[i]) / Cc[i] * (sigma22 - sigma11)/(np.log10(sigma22/sigma11))
                i +=1

        Me[0] = Me[1] /2 #adjust the most upper Me, because it must not be 0

        cv = self.ss.get_k() * Me / yw

        #length of vectors
        rows = len(self.ss.get_dz())
        #up = i-1; zero = i; lo = i+1
        up, zero, lo = self.ss.get_slices()
        #initialize
        fv = np.zeros((rows,))
        f1 = np.zeros((rows,))
        f2 = np.zeros((rows,))
        #factor vectors according to formula s.66
        fv[zero] = np.array(((1 + k[up] / k[zero]) / (1 + cv[zero] * k[up] / cv[up] / k[zero])) * cv[zero] * self.ss.dt / dz[up] ** 2)
        f1[zero] = np.array(2 * k[up] / (k[zero] + k[up]))
        f2[zero] = np.array(2 * k[zero] / (k[zero] + k[up]))
        #complete first and last entry
        fv[0], fv[-1] = fv[1], fv[-2]
        f1[0], f1[-1] = f1[1], f1[-2]
        f2[0], f2[-1] = f2[1], f2[-2]

        return fv, f1, f2


#solves drained-drained
    def solve00(self):

        #get the fixed data
        rows, A, up, zero, lo, cols, plottimes, plotmatrix, timelegend, ttrack, i = self.get_fixeddata()
        # get the factors
        fv, f1, f2 = self.get_factorvectors()
        #get drainvector
        drainvect = self.ss.get_drainvect()
        #get B and udisstot
        B, udisstot = A.copy(), A.copy()

        # FOR LOOP,
        for j in range(0, cols):

            #add load at time t
            for l, (time, load) in enumerate(self.tl):
                if ttrack == time:
                    A[zero] += load
            #internal drainage
            for j in drainvect:
                A[j] = self.dp #für drainage innerhalb der schichten A[] = 0 hier einfügen
            #save relevant vectors in plot matrix
            #damit alle dt funktionieren: add 'if i<len(timelegend)''
            if ttrack >= plottimes[i]:
                plotmatrix[:, [i]] = np.reshape(A,(rows,1))
                timelegend[[i]] = ttrack
                i += 1

            B = A.copy()
            # iteration zeitvektoren: CALCULATING NEXT TIME STEP
            # Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
            A[zero] = fv[zero] * (f1[zero] * A[up] - 2 * A[zero] + f2[zero] * A[lo]) + A[zero]

            deltau = B - A #ppressure change
            udisstot += deltau   # summed ppressure dissipated = summed sigmaeff change

            fv, f1, f2 = self.get_updated_factorvectors(deltau, udisstot)

            # timetracker: tt hat immer die Einheit der aktuellen Zeit in der Iteration (-> brauchbar für Zeiten des plots, und variable Lasten)
            ttrack += self.tt.dt
            ttrack = round(ttrack, 3)   #dies macht dt robuster (eliminates numerical noise which causes problems)

        return plotmatrix, timelegend, udisstot
#solves drained-undrained
    def solve01(self):

        #get the fixed data
        rows, A, up, zero, lo, cols, plottimes, plotmatrix, timelegend, ttrack, i = self.get_fixeddata()
        drainvect = self.ss.get_drainvect()

        # FOR LOOP,
        for j in range(0, cols):

            #get the variable data
            fv, f1, f2 = self.get_variabledata()

            #add load at time t
            for l, (time, load) in enumerate(self.tl):
                if ttrack == time:
                    A[zero] += load
                    A[-1] += load

            #internal drainage
            for j in drainvect:
                A[j] = self.dp #für drainage innerhalb der schichten A[] = 0 hier einfügen

            #save relevant vectors in plot matrix
                #damit alle dt funktionieren: add 'if i<len(timelegend)''
            if ttrack >= plottimes[i]:
                plotmatrix[:, [i]] = np.reshape(A,(rows,1))         #reshape entfernen?
                timelegend[[i]] = ttrack
                i += 1

            # iteration zeitvektoren:CALCULATING NEXT TIME STEP
            # Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
            A[zero] = fv[zero] * (f1[zero] * A[up] - 2 * A[zero] + f2[zero] * A[lo]) + A[zero]
            A[-1] = fv[-1] * (f1[-1] * A[-2] - 2 * A[-1] + f2[-1] * A[-2]) + A[-1]

            # timetracker: tt hat immer die Einheit der aktuellen Zeit in der Iteration (-> brauchbar für Zeiten des plots, und variable Lasten)
            ttrack += self.tt.dt
            ttrack = round(ttrack, 3)   #dies macht dt robuster (eliminates numerical noise which causes problems)

        return plotmatrix, timelegend
#solves undrained-drained
    def solve10(self):

        #get the fixed data
        rows, A, up, zero, lo, cols, plottimes, plotmatrix, timelegend, ttrack, i = self.get_fixeddata()
        drainvect = self.ss.get_drainvect()

        # FOR LOOP,
        for j in range(0, cols):

            #get the variable data
            fv, f1, f2 = self.get_variabledata()

            #add load at time t
            for l, (time, load) in enumerate(self.tl):
                if ttrack == time:
                    A[zero] += load
                    A[0] += load

            #internal drainage
            for j in drainvect:
                A[j] = self.dp #für drainage innerhalb der schichten A[] = 0 hier einfügen

            #save relevant vectors in plot matrix
                #damit alle dt funktionieren: add 'if i<len(timelegend)''
            if ttrack >= plottimes[i]:
                plotmatrix[:, [i]] = np.reshape(A,(rows,1))
                timelegend[[i]] = ttrack
                i += 1

            # iteration zeitvektoren:CALCULATING NEXT TIME STEP
            # Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
            A[zero] = fv[zero] * (f1[zero] * A[up] - 2 * A[zero] + f2[zero] * A[lo]) + A[zero]

            A[0] =   fv[0] * (f1[0] * A[1] - 2 * A[0] + f2[0] * A[1]) + A[0]

            # timetracker: tt hat immer die Einheit der aktuellen Zeit in der Iteration (-> brauchbar für Zeiten des plots, und variable Lasten)
            ttrack += self.tt.dt
            ttrack = round(ttrack, 3)   #dies macht dt robuster (eliminates numerical noise which causes problems)

        return plotmatrix, timelegend
#solves undrained-undrained
    def solve11(self):

        # get the fixed data
        rows, A, up, zero, lo, cols, plottimes, plotmatrix, timelegend, ttrack, i = self.get_fixeddata()
        drainvect = self.ss.get_drainvect()

        # FOR LOOP,
        for j in range(0, cols):

            # get the variable data
            fv, f1, f2 = self.get_variabledata()

            # add load at time t
            for l, (time, load) in enumerate(self.tl):
                if ttrack == time:
                    A[zero] += load
                    A[0] += load
                    A[-1] += load
            # internal drainage
            for j in drainvect:
                A[j] = self.dp  # für drainage innerhalb der schichten A[] = 0 hier einfügen

            # save relevant vectors in plot matrix
            # damit alle dt funktionieren: add 'if i<len(timelegend)''
            if ttrack >= plottimes[i]:
                plotmatrix[:, [i]] = np.reshape(A, (rows, 1))
                timelegend[[i]] = ttrack
                i += 1

            # iteration zeitvektoren:CALCULATING NEXT TIME STEP
            # Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
            A[zero] = fv[zero] * (f1[zero] * A[up] - 2 * A[zero] + f2[zero] * A[lo]) + A[zero]

            A[0] = fv[0] * (f1[0] * A[1] - 2 * A[0] + f2[0] * A[1]) + A[0]
            A[-1] = fv[-1] * (f1[-1] * A[-2] - 2 * A[-1] + f2[-1] * A[-2]) + A[-1]

            # timetracker: tt hat immer die Einheit der aktuellen Zeit in der Iteration (-> brauchbar für Zeiten des plots, und variable Lasten)
            ttrack += self.tt.dt
            ttrack = round(ttrack, 3)  # dies macht dt robuster (eliminates numerical noise which causes problems)

        return plotmatrix, timelegend

    def get_plot(self):

        # check mfactors !<0.5
        print('factors for each layer !<0.5')
        mfact = self.ss.get_mfact()
        print(mfact)
        #assert all(mfact < 0.5), 'check mfact, mathematically unstable'

        #solve while using the correct boundary conditions
        if self.bcs == [0, 0]:
            plotmatrix, timelegend, udisstot = self.solve00()
        elif self.bcs == [0, 1]:
            plotmatrix, timelegend = self.solve01()
        elif self.bcs == [1, 0]:
            plotmatrix, timelegend = self.solve10()
        elif self.bcs == [1, 1]:
            plotmatrix, timelegend = self.solve11()

        # plot erstellen
        #dzsum helps plotting without distortion
        dzsum = self.ss.get_dzsum()
        plt.plot(plotmatrix[:], -dzsum, label=timelegend)
        plt.xlabel("Pore water pressure")
        plt.ylabel("depth")
        plt.legend(loc=4, prop={'size': 6})
        plt.show()

        #print(udisstot)