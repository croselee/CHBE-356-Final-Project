import control
import numpy as np
from control import tf as tf
import control.matlab
import time

import glob
import pandas as pd

def read_profiles():
    filenames = glob.glob("Disturbance_profiles/profile*.csv")
    dfs = []
    for filename in filenames:
        df = pd.read_csv(filename, index_col=None, header=0)
        df = np.array(df)
        df = np.concatenate(df)
        dfs.append(df)
    return dfs

dfs = read_profiles()

tinitial = 0
tfinal = 59
N = len(dfs[0])
print(N)
times = np.linspace(tinitial, tfinal, N)
s = tf([1, 0], [0, 1])
Gp = 1 / (s ** 2 + s + 1)
Gd = (s + 1) / (s ** 2 + s + 1)

def calc_average_perf(s, Gp, Gd, Kc, tauI, tauD, tauC, times):
    Perf = np.linspace(0, 1, 100)
    Gc = Kc * (1 + 1 / (tauI * s) + tauD * s / (tauC * s + 1))
    sys_D = Gd / (1 + Gp * Gc)
    average_Perf = []
    for j in range(0, 100):
        Ti = dfs[j]
        _, Y, _ = control.forced_response(sys_D, times, Ti)
        sys_U = Gc
        _, U, _ = control.forced_response(sys_U, times, -Y)
        Perf[j] = sum(abs(Y)) + 0.2 * sum(abs(U))
        if sum(abs(Y) >= 5) >= 1 or sum(abs(U) >= 5) >= 1:
            Perf[j] = 1e6
        average_Perf = sum(Perf) / len(Perf)
    return average_Perf


def find_lowest_avg_perf():
    Npoints = 1
    Kc_start = 5.00E+03
    Kc_end = 9000
    tauI_start = 7.50E-01
    tauI_end = 9000
    tauD_start = 1.25E+03
    tauD_end = 9000
    tauC_start = 3.75E+03
    tauC_end = 9000

    avg_Perf = []
    for Kc in np.linspace(Kc_start, Kc_end, Npoints):
        for tauI in np.linspace(tauI_start, tauI_end, Npoints):
            for tauD in np.linspace(tauD_start, tauD_end, Npoints):
                for tauC in np.linspace(tauC_start, tauC_end, Npoints):
                    average_Perf = calc_average_perf(s, Gp, Gd, Kc, tauI, tauD, tauC, times)
                    print(Kc, tauI, tauD, tauC, average_Perf)
                    avg_Perf.append((Kc, tauI, tauD, tauC, average_Perf))


from scipy.optimize import minimize

def find_lowest_avg_perf_2():
    x0 = [5.00E+03, 7.50E-01, 1.25E+03, 3.75E+03]

    def f(params):
        Kc, tauI, tauD, tauC = params
        return calc_average_perf(s, Gp, Gd, Kc, tauI, tauD, tauC, times)

    print(minimize(f, x0))

def only_one_MV():
    Npoints = 100
    Kc = 5
    tauI = 5
    tauD = 5
    tauC = 5
    MV_start = 1
    MV_end = 5000

    avg_Perf = []
    for manipulated in np.linspace(MV_start, MV_end, Npoints):
        tauC = manipulated
        average_Perf = calc_average_perf(s, Gp, Gd, Kc, tauI, tauD, tauC, times)
        print(Kc, tauI, tauD, tauC, average_Perf)
        avg_Perf.append((Kc, tauI, tauD, tauC, average_Perf))

def only_two_MV():
    Npoints = 20
    Kc_start = 1
    Kc_end = 6000
    tauI = 100
    tauD_start = 1
    tauD_end = 6000
    tauC = 100

    avg_Perf = []
    for Kc in np.linspace(Kc_start, Kc_end, Npoints):
        for tauD in np.linspace(tauD_start, tauD_end, Npoints):
               average_Perf = calc_average_perf(s, Gp, Gd, Kc, tauI, tauD, tauC, times)
               print(Kc, tauI, tauD, tauC, average_Perf)
               avg_Perf.append((Kc, tauI, tauD, tauC, average_Perf))

now = time.time()
#find_lowest_avg_perf()
find_lowest_avg_perf_2()
#only_one_MV()
#only_two_MV()
print(now - time.time())

