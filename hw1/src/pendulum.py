import numpy as np
import matplotlib.pyplot as plt
import math
import os

"""
The equation we want to solve:
\theta(t+h)=2\theta(t)-\theta(t-h)+(-g/l)\theta(t)h^2
\theta(t=0)=1/30
"""


def mkdir(path):

    folder = os.path.exists(path)

    if not folder:
        os.makedirs(path)
        print("create new folder")
    else:
        print("folder exists")


def theta_generator(theta_0, h, time_length):
    # initial condition
    i, theta_t_minus_h, theta_t = 0, theta_0, theta_0
    while i < time_length:
        yield theta_t
        # Verlet method
        if theta_0 < float(10 / 180 * math.pi):
            theta_t_plus_h = 2 * theta_t - theta_t_minus_h - g / L * theta_t * h**2
        else:
            theta_t_plus_h = 2 * theta_t - theta_t_minus_h - g / L * math.sin(theta_t) * h**2
        # propagte
        i, theta_t, theta_t_minus_h = i + 1, theta_t_plus_h, theta_t


def compute_theta_energy(theta_0, h, time_length):
    # theta list with time
    theta = [x for x in theta_generator(theta_0, h, time_length)]
    # the unit of time is `h` s
    time = [t for t in range(time_length)]

    # theta + h/2 with time, half represent time at t + h/2
    time_half = [t + h / 2 for t in time]
    time_half.pop(-1)
    theta_half = [(theta[i + 1] + theta[i]) / 2 for i in range(len(theta) - 1)]
    velocity_half = [L * (theta[i + 1] - theta[i]) / (2 * h) for i in range(len(theta) - 1)]
    T_half = [1 / 2 * m * v**2 for v in velocity_half]
    V_half = [m * g * L * (1 - math.cos(theta)) for theta in theta_half]
    E_half = np.array(T_half) + np.array(V_half)
    assert len(T_half) == len(V_half), "Error"
    return time, theta, time_half, E_half


def make_plots(time, theta, time_half, E_half, E_exact):
    mkdir("./figs")
    plt.plot(time, theta)
    plt.savefig("./figs/theta_time.png", dpi=400)
    plt.clf()

    plt.plot(time_half, E_half)
    plt.hlines(E_exact, time_half[0], time_half[-1], colors="c", linestyles="dashed")
    plt.savefig("./figs/energy_time.png", dpi=400)


if __name__ == "__main__":
    # constant
    g = 9.8
    L = 0.1
    h = 0.01
    m = 1
    time_length = 100

    # initial condition
    theta_0 = float(1 / 30)
    # exact energy
    E_exact = m * g * L * (1 - math.cos(theta_0))
    time, theta, time_half, E_half = compute_theta_energy(theta_0, h, time_length)
    make_plots(time, theta, time_half, E_half, E_exact)
