
import timee as tm
import assembly as am
import layer as lm
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

class Model:
    def __init__(self, tl, ss, tt, graphs, dp, yw):
        self.tl = tl    #ja
        self.ss = ss    #ja
        self.tt = tt    #ja
        self.graphs = graphs
        self.dp = dp
        self.yw = yw

    def get_fixeddata(self):    #this prepares the fix data for the iteration

        # check mfactors !<0.5
        print('factors for each layer !<0.5')
        mfact = self.ss.get_mfact()
        print(mfact)
        assert all(mfact < 0.5), 'check mfact, mathematically unstable'

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

    def get_factor_fun(self):
        #needed vectors
        Me = self.ss.get_Me()
        dz = self.ss.get_dz()
        k = self.ss.get_k()
        e0 = self.ss.get_e0()
        Cc = self.ss.get_Cc()
        dt = self.ss.dt
        effsigma0 = self.ss.get_effsigma()

        rows = len(self.ss.get_dz())
        # up = i-1; zero = i; lo = i+1
        up, zero, lo = self.ss.get_slices()

        def fun(deltau, udisstot):
            # calculate sigma1 (vor it.schritt) and sigma2 (nach it.schritt)
            sigma2 = effsigma0 + udisstot
            sigma1 = sigma2 - deltau

            # calculate new Me
            Me = np.log(10) * (1 + e0) / Cc * sigma2

            Me[0] = Me[1] / 2  # adjust the most upper Me, because it must not be 0

            cv = k * Me / self.yw

            # initializen factor vectors
            fv, f1, f2 = np.zeros((rows,)), np.zeros((rows,)), np.zeros((rows,))
            # factor vectors according to formula s.66
            fv[zero] = np.array(((1 + k[up] / k[zero]) / (1 + cv[zero] * k[up] / cv[up] / k[zero])) * cv[zero] * dt / dz[up] ** 2)
            f1[zero] = np.array(2 * k[up] / (k[zero] + k[up]))
            f2[zero] = np.array(2 * k[zero] / (k[zero] + k[up]))
            # complete first and last entry
            fv[0], fv[-1] = fv[1], fv[-2]
            f1[0], f1[-1] = f1[1], f1[-2]
            f2[0], f2[-1] = f2[1], f2[-2]
            return fv, f1, f2
        return fun


#solves all
    def solve(self, top_drained=True, bot_drained=True, non_linear=True):

        #get the fixed data
        rows, A, up, zero, lo, cols, plottimes, plotmatrix, timelegend, ttrack, i = self.get_fixeddata()
        # get the factors
        fv, f1, f2 = self.get_factorvectors()
        #get drainvector
        drainvect = self.ss.get_drainvect()
        factor_fun = self.get_factor_fun()

        #get B and udisstot
        B, udisstot = A.copy(), A.copy()

        mask_load = slice(top_drained, rows-bot_drained)
        # FOR LOOP,
        for j in range(0, cols):

            #add load at time t
            for time, load in self.tl:
                if ttrack == time:
                    A[mask_load] += load

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

            if not bot_drained:
                A[-1] = fv[-1] * (f1[-1] * A[-2] - 2 * A[-1] + f2[-1] * A[-2]) + A[-1]
            if not top_drained:
                A[0] = fv[0] * (f1[0] * A[1] - 2 * A[0] + f2[0] * A[1]) + A[0]

            #internal drainage do it again, so the porewaterpressure doesnt increase (very bad for nonlinear solution)
            for j in drainvect:
                A[j] = self.dp #für drainage innerhalb der schichten A[] = 0 hier einfügen

            deltau = B - A #ppressure change
            udisstot += deltau   # summed ppressure dissipated = summed sigmaeff change

            fv, f1, f2 = factor_fun(deltau, udisstot)

            # timetracker: tt hat immer die Einheit der aktuellen Zeit in der Iteration (-> brauchbar für Zeiten des plots, und variable Lasten)
            ttrack += self.tt.dt
            ttrack = round(ttrack, 3)   #dies macht dt robuster (eliminates numerical noise which causes problems)

        return Solution(self.ss, plottimes, plotmatrix)


class Solution:
    def __init__(self, assembly, times, pore_pressures):
        self.assembly = assembly
        self.times = times
        self.pore_pressures = pore_pressures

    def plot_pressures(self, plot_times, plot_depths=None):
        if plot_depths is None:
            plot_depths = self.assembly.get_dzsum()

        plot_pressures = interp2d(self.times, self.assembly.get_dzsum(), self.pore_pressures)(plot_times, plot_depths)
        plt.plot(plot_pressures, -plot_depths, label=plot_times)
        plt.xlabel("Pore water pressure")
        plt.ylabel("depth")
        plt.title('u(t)-z-diagram')
        plt.legend(loc=4, prop={'size': 6})
        plt.show()
        return

    def plot_U(self):

        ut0 = sum(self.pore_pressures[: , 0])
        Uvector = self.times.copy()
        wi = self.assembly.get_dz() #weil dz unterschiedlich ist, muss gewichtet werden bzgl dz
        avgwi = np.mean(wi, (0,))

        i = 0
        for counter in self.times:
            ut = sum(self.pore_pressures[:, i] * wi/avgwi)
            Uvector[i] = ut/ut0
            i += 1

        plt.plot(self.times, 1-Uvector, label = 'U')
        plt.ylim([0, 1])
        plt.gca().invert_yaxis()
        plt.xlabel("Time [s]")
        plt.ylabel("avg. consolidation U")
        plt.title('U(t)')
        plt.legend(loc=1, prop={'size': 6})
        plt.show()

#settlement simple approach (top and bottom layer are not calculated)
    def plot_settlement(self):

        um =       self.pore_pressures
        uincrm =   np.zeros(um.shape)
        sigeffm =  np.zeros(um.shape)
        evolm =    np.zeros(um.shape)
        settlemm = np.zeros(um.shape)
        sigeff0 =  self.assembly.get_effsigma()
        settlemvect = np.zeros(self.times.shape)

        i = 1
        for counter in um[0, 1:]: #calculate uincrements
            if any(um[:, i-1] - um[:, i] < 0):
                uincrm[:, i] = 0
            else:
                uincrm[:, i] = um[:, i-1] - um[:, i]
            i += 1

        sigeffm[:, 0] = sigeff0
        i = 1
        for counter in um[0, 1:]:   #calculate sigeff at time t
            sigeffm[:, i] = sigeffm[:, i-1] + uincrm[:, i]
            i += 1

        i=0
        Cc = self.assembly.get_Cc()
        e = self.assembly.get_e0()
        for counter in um[0,:]:     #calculate strains
            evolm[1:, i] = Cc[1:] * np.log10(sigeffm[1:, i]/sigeffm[1:, 0]) / (1 + e[1:])
            i += 1

        i=0
        for counter in um[0, :]:    #calculate settlement = dz*epsilon
            settlemm[:, i] = evolm[:, i] * self.assembly.get_dz()
            i+=1

        ##korrektur: first and last layer linear erweitern mittels 2&3 eintrag
        #settlemm[0,:] = settlemm[1,:] + settlemm[1,:] - settlemm[2,:]
        #settlemm[-1,:] = settlemm[-2,:] + settlemm[-2,:] - settlemm[-3,:]

        i=0
        for counter in um[0, :]:
            settlemvect[i] = sum(settlemm[:, i])
            i += 1

        plt.plot(self.times, -settlemvect, label = 'Settlement [m]')
        plt.xlabel("Time [s]")
        plt.ylabel("Settlement[m]")
        plt.title('Settlement(t) 1st version')
        plt.legend(loc=1, prop={'size': 6})
        plt.show()

        return sigeffm, evolm

#settlement interpolated approach
    def plot_settlement2(self, tl, top_drained=True, bot_drained=True):
        self.tl = tl

        sigeffm, evolm = self.plot_settlement()
        sigeffmavg = np.zeros(sigeffm.shape)

        #immediately add the load at time t in drained points
        i=0
        for time in self.times:
            for t, l in tl:
                if t <= time:
                    if top_drained:
                        sigeffm[0, i] += l
                    if bot_drained:
                        sigeffm[-1, i] += l
            i+=1

        i=0
        for counter in sigeffmavg[:-1,0]:
            sigeffmavg[i,:] = (sigeffm[i,:] + sigeffm[i+1,:]) /2
            i+=1

        sigeff0 = self.assembly.get_effsigma()
        sigeff0avg = np.zeros(sigeff0.shape)

        i=0
        for counter in sigeff0[:-1]:
            sigeff0avg[i] = (sigeff0[i] + sigeff0[i+1]) /2
            i+=1

        evolmavg = np.zeros(sigeffmavg.shape)

        Cc = self.assembly.get_Cc()
        e = self.assembly.get_e0()

        i = 0
        for counter in evolmavg[0,:]:     #nur bis :-1 weil der letzte eintrag nicht nötig ist.
            evolmavg[:-1, i] = Cc[:-1] * np.log10(sigeffmavg[:-1, i]/sigeff0avg[:-1]) / (1 + e[:-1])
            i += 1

        settlemmavg = np.zeros(evolm.shape)
        i=0
        for counter in settlemmavg[0, :]:    #calculate settlement = dz*epsilon
            settlemmavg[:, i] = evolmavg[:, i] * self.assembly.get_dz()
            i+=1

        settlemvectavg = np.zeros(self.times.shape)

        i=0
        for counter in settlemmavg[0, :]:
            settlemvectavg[i] = sum(settlemmavg[:, i])
            i += 1

        plt.plot(self.times, -settlemvectavg, label = 'Settlement [m]')
        plt.xlabel("Time [s]")
        plt.ylabel("Settlement[m]")
        plt.title('Settlement(t) 2nd version')
        plt.legend(loc=1, prop={'size': 6})
        plt.show()

    def get_plot(self, **kwargs):

        # check mfactors !<0.5
        print('factors for each layer !<0.5')
        mfact = self.assembly.get_mfact()
        print(mfact)
        #assert all(mfact < 0.5), 'check mfact, mathematically unstable'

        #solve while using the correct boundary conditions
        plotmatrix, timelegend = self.pore_pressures, self.times

        # plot erstellen
        #dzsum helps plotting without distortion
        dzsum = self.assembly.get_dzsum()
        plt.plot(plotmatrix[:], -dzsum, label=timelegend)
        plt.xlabel("Pore water pressure")
        plt.ylabel("depth")
        plt.legend(loc=4, prop={'size': 6})
        plt.show()

        #print(udisstot)
