import control
import numpy as np
from scipy.interpolate import interp1d

NUM_PROFILES = 100
MAX_TIME = 60
T_sp = 0
K_c, tau_I, tau_D, tau_c = [6934.46250453,  110.42637161, 4551.44114479, 9952.96261485]

def simulate(i, K_c, tau_I, tau_D, tau_c):
    time_data = np.linspace(0, MAX_TIME-1, MAX_TIME)
    with open(f'profile{i:02}.csv', 'r') as f:
        f.readline()    # skip header row
        T_i_data = [float(line.strip()) for line in f]
    f = interp1d(time_data, T_i_data)
    points_per_seond = 100
    time = np.linspace(0, MAX_TIME-1, MAX_TIME*points_per_seond)
    T_i = f(time)

    s = control.tf([1,0],[0,1])
    G_p = 1/(s**2 + s + 1)
    G_d = (s+1)/(s**2 + s + 1)
    G_c = K_c * (1+ 1/(tau_I*s) + (tau_D*s)/(tau_c*s + 1))
    sys_D = G_d / (1 + G_p * G_c)

    Tsp = 0
    _, T, _ = control.forced_response(sys_D, time, T_i)
    _, Q, _ = control.forced_response(G_c, time, Tsp - T)
    error = (sum(abs(T)) + 0.2*sum(abs(Q)))/points_per_seond
    if sum(abs(T)>= 5) or sum(abs(Q)>= 5):
        error = 1e6
    return error

total_error = 0
for i in range(NUM_PROFILES):
    error = simulate(i, K_c, tau_I, tau_D, tau_c)
    total_error += error
print(f'Average error: {total_error / NUM_PROFILES}')