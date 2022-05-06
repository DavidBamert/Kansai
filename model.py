
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

yw = 10


class Model:
    def __init__(self, tl, ss, tt, graphs, dp):
        self.tl = tl
        self.ss = ss
        self.tt = tt
        self.graphs = graphs
        self.dp = dp

    # This function prepares the fix data for the iteration
    def get_fixeddata(self):

        """
        # Check mfactors !<0.5
        print('factors for each layer !<0.5')
        mfact = self.ss.get_mfact()
        print(mfact)
        assert all(mfact < 0.5), 'check mfact, mathematically unstable'
        """

        # Create A
        rows = len(self.ss.get_dz())
        A = np.zeros((rows,))

        # Get slices
        up, zero, lo = self.ss.get_slices()

        # Get time discretization and plot
        cols = self.tt.get_cols()
        plottimes, plotmatrix, timelegend = self.tt.get_plotmatrix(rows, self.graphs)  # ev. without graphs

        # ttrack = time tracker
        ttrack = 0
        # i = column number of plotmatrix
        i = 0

        return rows, A, up, zero, lo, cols, plottimes, plotmatrix, timelegend, ttrack, i

    # This function prepares (and in the future alters) the variable data for the iteration
    def get_factorvectors(self):
        # Change fi and dz
        fv, f1, f2 = self.ss.get_factors()
        return fv, f1, f2

    def get_updated_factorvectors(self, deltau, udisstot):
        # Needed vectors
        Me = self.ss.get_Me()
        dz = self.ss.get_dz()
        k = self.ss.get_k()
        e0 = self.ss.get_e0()
        Cc = self.ss.get_Cc()
        # Calculate sigma1 and sigma2
        sigma2 = self.ss.get_effsigma() + udisstot
        sigma1 = sigma2 - deltau
        # Calculate new e
        j = 0
        for sigma11, sigma22 in zip(sigma1, sigma2):  # wie kann man hier fallunterscheiden? !!delta e0 aus jeder iteration müsste summiert sein!!
            if sigma22-sigma11 > 1e-7:
                e0[j] = e0[j] - Cc[j] * np.log10(sigma22/sigma11)
                j += 1

        # Calculate new Me
        Me = np.log(10) * (1 + e0) / Cc * sigma2
        """
        #tangentenmodul ist immer besser (?)
        i = 0
        for sigma11,sigma22 in zip(sigma1,sigma2): #problem with log <- fallunterscheidung für die Me berechnung
            if sigma22-sigma11 < 1e-7:
                Me[i] = np.log(10) * (1+e0[i]) / Cc[i] * sigma22
                i +=1
            else:
                Me[i] = (1 + e0[i]) / Cc[i] * (sigma22 - sigma11)/(np.log10(sigma22/sigma11))
                i +=1
        """

        Me[0] = Me[1]/2  # adjust the most upper Me, because it must not be 0

        cv = self.ss.get_k() * Me / yw

        # Length of vectors
        rows = len(self.ss.get_dz())
        # up = i-1; zero = i; lo = i+1
        up, zero, lo = self.ss.get_slices()
        # Initialize
        fv, f1, f2 = np.zeros((rows,)), np.zeros((rows,)), np.zeros((rows,))
        # Factor vectors according to formula s.66
        fv[zero] = np.array(((1 + k[up] / k[zero]) / (1 + cv[zero] * k[up] / cv[up] / k[zero])) * cv[zero] * self.ss.dt / dz[up] ** 2)
        f1[zero] = np.array(2 * k[up] / (k[zero] + k[up]))
        f2[zero] = np.array(2 * k[zero] / (k[zero] + k[up]))
        # Complete first and last entry
        fv[0], fv[-1] = fv[1], fv[-2]
        f1[0], f1[-1] = f1[1], f1[-2]
        f2[0], f2[-1] = f2[1], f2[-2]

        return fv, f1, f2

    # Solves all
    def solve(self, top_drained=True, bot_drained=True):

        # Get the fixed data
        rows, A, up, zero, lo, cols, plottimes, plotmatrix, timelegend, ttrack, i = self.get_fixeddata()
        # Get the factors
        fv, f1, f2 = self.get_factorvectors()
        # Get drainvector
        drainvect = self.ss.get_drainvect()
        # Get B and udisstot
        B, udisstot = A.copy(), A.copy()

        mask_load = slice(top_drained, rows-bot_drained)  # False = 0, True = 1
        # FOR LOOP,
        for j in range(0, cols):

            # Add load at time t
            for time, load in self.tl:
                if ttrack == time:
                    A[mask_load] += load

            # Internal drainage
            for m in drainvect:
                A[m] = self.dp  # für drainage innerhalb der schichten A[] = 0 hier einfügen

            # Save relevant vectors in plot matrix
            # Damit alle dt funktionieren: add 'if i<len(timelegend)''
            if ttrack >= plottimes[i]:
                plotmatrix[:, [i]] = np.reshape(A, (rows, 1))
                timelegend[[i]] = ttrack
                i += 1

            B = A.copy()
            # Iteration zeitvektoren: CALCULATING NEXT TIME STEP
            # Übergangsbedingung eignet sich als allgemeinere Formel! (Buch s.66)
            A[zero] = fv[zero] * (f1[zero] * A[up] - 2 * A[zero] + f2[zero] * A[lo]) + A[zero]

            if not bot_drained:
                A[-1] = fv[-1] * (f1[-1] * A[-2] - 2 * A[-1] + f2[-1] * A[-2]) + A[-1]
            if not top_drained:
                A[0] = fv[0] * (f1[0] * A[1] - 2 * A[0] + f2[0] * A[1]) + A[0]

            deltau = B - A  # ppressure change
            udisstot += deltau  # summed ppressure dissipated = summed sigmaeff change

            fv, f1, f2 = self.get_updated_factorvectors(deltau, udisstot)

            # timetracker: tt hat immer die Einheit der aktuellen Zeit in der Iteration (-> brauchbar für Zeiten des plots, und variable Lasten)
            ttrack += self.tt.dt
            ttrack = round(ttrack, 3)   # dies macht dt robuster (eliminates numerical noise which causes problems)

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
        plt.xlabel("Excess pore water pressure [kPa]")
        plt.ylabel("Depth [m]")
        plt.title('u(t)-z-diagram')
        plt.legend(loc=4, prop={'size': 6})
        plt.show()
        return

    def get_U(self):

        ut0 = sum(self.pore_pressures[:, 0])
        Uvector = self.times.copy()
        wi = self.assembly.get_dz()  # weil dz unterschiedlich ist, muss gewichtet werden bzgl dz
        avgwi = np.mean(wi, (0,))

        i = 0
        for counter in self.times:
            ut = sum(self.pore_pressures[:, i] * wi/avgwi)
            Uvector[i] = ut/ut0
            i += 1

        plt.plot(self.times, Uvector, label='U')
        plt.xlabel("Time [s]")
        plt.ylabel("avg. consolidation U")
        plt.title('U(t)')
        plt.legend(loc=1, prop={'size': 6})
        plt.show()

    def get_settlement(self):

        effsigma, effsigma0 = self.assembly.get_effsigma(), self.assembly.get_effsigma()

        ut1 = sum(self.pore_pressures[i , 0])
        ut2 = sum(self.pore_pressures[i+1 , 0])
        uincr = ut1-ut2

        if any(uincr) < 0:
            uincr[i] = 0
        else:
            uincr

        strain = self.assembly.get_Cc()
        de = -self.assembly.get_Cc() * np.log10(effsigma+uincr/effsigma)

        return 0

    def get_dzz(self):
        z = sum(self.assembly.get_dz()[0:-1])
        print(z)
        return z


""" ancient relic
    def get_plot(self, **kwargs):

        # check mfactors !<0.5
        print('factors for each layer !<0.5')
        mfact = self.ss.get_mfact()
        print(mfact)
        #assert all(mfact < 0.5), 'check mfact, mathematically unstable'

        #solve while using the correct boundary conditions
        plotmatrix, timelegend, udisstot = self.solve(**kwargs)

        # plot erstellen
        #dzsum helps plotting without distortion
        dzsum = self.ss.get_dzsum()
        plt.plot(plotmatrix[:], -dzsum, label=timelegend)
        plt.xlabel("Pore water pressure")
        plt.ylabel("depth")
        plt.legend(loc=4, prop={'size': 6})
        plt.show()

        #print(udisstot)
"""