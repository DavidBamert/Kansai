import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d


class Model:
    def __init__(self, tl, ss, tt, graphs, dp, yw):
        self.tl = tl
        self.ss = ss
        self.tt = tt
        self.graphs = graphs
        self.dp = dp
        self.yw = yw

    def get_fixeddata(self):

        # create A
        rows = len(self.ss.get_dz())
        A = np.zeros((rows,))
        # get factors and slices <-bcs
        up, zero, lo = self.ss.get_slices()
        # get time discretization and plot
        cols = self.tt.get_cols()
        plottimes, plotmatrix, timelegend = self.tt.get_plotmatrix(rows, self.graphs)
        dt = self.tt.dt

        return rows, A, up, zero, lo, cols, plottimes, plotmatrix, timelegend, dt

    # function for non-linear calculation
    def get_factor_fun(self):

        # needed vectors
        rows = len(self.ss.get_dz())
        dz = self.ss.get_dz()
        k = self.ss.get_k()
        e0 = self.ss.get_e0()
        Cc = self.ss.get_Cc()
        dt = self.ss.dt
        effsigma0 = self.ss.get_effsigma0()
        # initializen factor vectors
        fv, f1, f2 = np.zeros((rows,)), np.zeros((rows,)), np.zeros((rows,))
        # up = i-1; zero = i; lo = i+1
        up, zero, lo = self.ss.get_slices()

        # 2nd order dz
        effsigma0avg = np.zeros(effsigma0.shape)
        effsigmatavg = np.zeros(effsigma0.shape)
        epsvolavg = np.zeros(effsigma0.shape)
        epsvol = np.zeros(effsigma0.shape)

        def fun(udisstot, sec_order_strains):
            # calculate sigma1 (before it.step) and sigma2 (after it.step)
            sigma2 = effsigma0 + udisstot

            # calculate new Me
            Me = np.log(10) * (1 + e0) / Cc * sigma2
            Me[0] = Me[1] / 2  # adjust the most upper Me, because it must not be 0

            cv = k * Me / self.yw
            dzsecnd = dz

            if sec_order_strains:
                # calculate new dz; last entry always 0, because useless
                effsigma0avg[0:-1] = (effsigma0[0:-1] + effsigma0[1:]) / 2
                effsigmatavg[0:-1] = (sigma2[0:-1] + sigma2[1:]) / 2
                epsvolavg[0:-1] = Cc[0:-1] * np.log10(effsigmatavg[0:-1] / effsigma0avg[0:-1]) / (1 + e0[0:-1])
                # return from averages to discretization points
                epsvol[1:-1] = (epsvolavg[0:-2] + epsvolavg[1:-1]) / 2
                # complete first and last entry <- assumption
                epsvol[0], epsvol[-1] = epsvolavg[0], epsvolavg[-1]
                dzsecnd = dz - (epsvol * dz)
                fv[zero] = np.array(((1 + k[up] / k[zero]) / (1 + cv[zero] * k[up] / cv[up] / k[zero])) * cv[zero] * dt / dzsecnd[up] ** 2)
            if not sec_order_strains:
                fv[zero] = np.array(((1 + k[up] / k[zero]) / (1 + cv[zero] * k[up] / cv[up] / k[zero])) * cv[zero] * dt / dz[up] ** 2)

            f1[zero] = np.array(2 * k[up] / (k[zero] + k[up]))
            f2[zero] = np.array(2 * k[zero] / (k[zero] + k[up]))
            # complete first and last entry
            fv[0], fv[-1] = fv[1], fv[-2]
            f1[0], f1[-1] = f1[1], f1[-2]
            f2[0], f2[-1] = f2[1], f2[-2]

            return fv, f1, f2, cv, dzsecnd
        return fun


    # solves all
    def solve(self, top_drained=True, bot_drained=True, non_linear=True, sec_order_strains=True):

        # check mfactors !<0.5 before START LOOP
        if non_linear:
            mfact0, maxdepth = self.ss.get_mfact0_nonlin()
            print('At depth', maxdepth, 'm, max. M factor:')
        else:
            mfact0 = self.ss.get_mfact0_lin()
            print('M factors for each layer')
        print(mfact0)
        assert all(mfact0 < 0.5), 'Initial mfact > 0.5, mathematically unstable!'

        # get the fixed data
        rows, A, up, zero, lo, cols, plottimes, plotmatrix, timelegend, dt = self.get_fixeddata()
        udisstot = A.copy()
        mask_load = slice(top_drained, rows - bot_drained)

        # get the initial factors for the first calculation of A
        fv, f1, f2 = self.ss.get_factors0(non_linear)

        factor_fun = self.get_factor_fun()
        drainvect = self.ss.get_drainvect()

        ttrack = 0  # ttrack = timetracker
        i = 0  # i = column number of plotmatrix

        # START LOOP
        for j in range(0, cols):

            # add load at time t
            for time, load in self.tl:
                if ttrack == time:
                    A[mask_load] += load
                    if bot_drained:  # also at time t, add the load as dissipated pore pressure at drained locations
                        udisstot[-1] += load
                    if top_drained:
                        udisstot[0] += load
                    for j in drainvect:
                        udisstot[j] += load

            # internal drainage
            for j in drainvect:
                A[j] = self.dp  # for drainage in layer A[] = 0 insert here

            # save relevant vectors in plot matrix
            # so that all dt works: add 'if i<len(timelegend)'
            if ttrack >= plottimes[i]:
                plotmatrix[:, [i]] = np.reshape(A, (rows, 1))
                timelegend[[i]] = ttrack
                i += 1

            B = A.copy()
            # iteration time vectors: CALCULATING NEXT TIME STEP
            # transition condition is suitable as a more general formula! (LHAP page 66)
            A[zero] = fv[zero] * (f1[zero] * A[up] - 2 * A[zero] + f2[zero] * A[lo]) + A[zero]
            if not bot_drained:
                A[-1] = fv[-1] * (f1[-1] * A[-2] - 2 * A[-1] + f2[-1] * A[-2]) + A[-1]
            if not top_drained:
                A[0] = fv[0] * (f1[0] * A[1] - 2 * A[0] + f2[0] * A[1]) + A[0]

            # internal drainage do it again, so the pore water pressure doesn't increase (very bad for nonlinear solution)
            for j in drainvect:
                A[j] = self.dp  # for drainage in layer A[] = 0 insert here

            deltau = B - A  # ppressure change
            udisstot += deltau   # summed ppressure dissipated = summed sigmaeff change

            if non_linear:
                fv, f1, f2, cv, dzsecnd = factor_fun(udisstot, sec_order_strains)

            # timetracker: tt always has the unit of the current time in the iteration (-> useful for times of plots, and variable loads)
            ttrack += dt
        # END LOOP

        # check mfactors !<0.5 after END LOOP (cv have increased in non-linear case)
        if non_linear:
            mfact = cv * dt / dzsecnd ** 2
            maxfact = max(mfact)
            for i, fact in enumerate(mfact):
                if fact == maxfact:
                    maxdepth = self.ss.get_dzsum()[i]
            assert all(mfact < 0.5), f"At depth {maxdepth} m, maximum final mfact {maxfact}!"

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
        plt.title('u(t)-z diagram')
        plt.legend(loc=4, prop={'size': 6})
        plt.show()
        return -plot_pressures

    def plot_U(self):

        ut0 = sum(self.pore_pressures[:, 0])
        Uvector = self.times.copy()
        wi = self.assembly.get_dz()  # because dz is different, it must be weighted relative to dz
        avgwi = np.mean(wi, (0,))

        i = 0
        for counter in self.times:
            ut = sum(self.pore_pressures[:, i] * wi/avgwi)
            Uvector[i] = ut/ut0
            i += 1

        plt.plot(self.times, 1-Uvector, label='U')
        plt.ylim([0, 1])
        plt.gca().invert_yaxis()
        plt.xlabel("Time [s]")
        plt.ylabel("Avg. consolidation U")
        plt.title('U(t)')
        plt.legend(loc=1, prop={'size': 6})
        plt.show()


    # settlement interpolated approach
    def plot_settlement(self, tl, top_drained=True, bot_drained=True):
        self.tl = tl

        um = self.pore_pressures
        uincrm = np.zeros(um.shape)
        sigeffm = np.zeros(um.shape)
        sigeff0 = self.assembly.get_effsigma0()
        sigeff0avg = np.zeros(sigeff0.shape)

        i = 0
        for counter in sigeff0[:-1]:
            sigeff0avg[i] = (sigeff0[i] + sigeff0[i + 1]) / 2
            i += 1

        i = 1
        for counter in um[0, 1:]:  # calculate uincrements
            if any(um[:, i - 1] - um[:, i] < 0):
                uincrm[:, i] = 0
            else:
                uincrm[:, i] = um[:, i - 1] - um[:, i]
            i += 1

        sigeffm[:, 0] = sigeff0
        i = 1
        for counter in um[0, 1:]:  # calculate sigeff at time t
            sigeffm[:, i] = sigeffm[:, i - 1] + uincrm[:, i]
            i += 1
        # immediately add the load at time t (and forward) in drained points
        i = 0
        for time in self.times:
            for t, l in tl:
                if t <= time:
                    if top_drained:
                        sigeffm[0, i] += l
                    if bot_drained:
                        sigeffm[-1, i] += l
                    for j in self.assembly.get_drainvect():
                        sigeffm[j, i] += l
            i += 1

        sigeffmavg = np.zeros(sigeffm.shape)
        i = 0
        for counter in sigeffmavg[:-1, 0]:
            sigeffmavg[i, :] = (sigeffm[i, :] + sigeffm[i+1, :]) / 2
            i += 1

        evolmavg = np.zeros(sigeffmavg.shape)
        i = 0
        Cc = self.assembly.get_Cc()
        e = self.assembly.get_e0()
        for counter in evolmavg[0, :]:  # only up to :-1 because the last entry is not necessary
            evolmavg[:-1, i] = Cc[:-1] * np.log10(sigeffmavg[:-1, i]/sigeff0avg[:-1]) / (1 + e[:-1])
            i += 1

        settlemmavg = np.zeros(um.shape)
        i = 0
        for counter in settlemmavg[0, :]:  # calculate settlement = dz*epsilon
            settlemmavg[:, i] = evolmavg[:, i] * self.assembly.get_dz()
            i += 1

        settlemvectavg = np.zeros(self.times.shape)

        i = 0
        for counter in settlemmavg[0, :]:
            settlemvectavg[i] = sum(settlemmavg[:, i])
            i += 1

        plt.plot(self.times, -settlemvectavg, label='Settlement [m]')
        plt.ylim([0, -settlemvectavg[-1] - 0.1])
        plt.gca().invert_yaxis()
        plt.xlabel("Time [s]")
        plt.ylabel("Settlement [m]")
        plt.title('Settlement(t)')
        plt.legend(loc=1, prop={'size': 6})
        plt.show()

        return -np.reshape(settlemvectavg, (len(settlemvectavg), 1))
