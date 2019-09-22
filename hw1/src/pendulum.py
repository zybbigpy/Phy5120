import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
import os


def theta_generator(theta_0, h, time_length):
    """
    The equation we want to solve:
    \theta(t+h)=2\theta(t)-\theta(t-h)+(-g/l)\theta(t)h^2
    \theta(t=0)=1/30
    """
    # initial condition
    i, theta_t_minus_h, theta_t = 0, theta_0, theta_0
    while i < time_length:
        yield theta_t
        # Verlet method
        if theta_0 < float(10 / 180 * math.pi):  # for small angle
            theta_t_plus_h = 2 * theta_t - theta_t_minus_h - g / L * theta_t * h**2
        else:  # for large angle
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


def make_plots(time, theta, time_half, E_half, E_exact, h, theta_0):

    plt.plot(time, theta)
    plt.ylabel("$\\theta$")
    plt.xlabel("Time: " + str(h) + "s")
    plt.title("$\\theta$-time relationship")
    plt.savefig("./figs/theta_time_" + str(h) + "_" + str(theta_0) + ".png", dpi=400)
    plt.clf()

    plt.plot(time_half, E_half)
    plt.xlabel("Time:" + str(h) + "s")
    plt.ylabel("Energy: J")
    plt.xlim(time[0], time[-1] + 50)
    plt.hlines(E_exact, time_half[0], time_half[-1], colors="c", linestyles="dashed")
    plt.text(time_half[-1], E_exact, 'Exact Energy ', ha='left', va='center')
    plt.title("Energy-time relationship")
    plt.savefig("./figs/energy_time_" + str(h) + "_" + str(theta_0) + ".png", dpi=400)


if __name__ == "__main__":
    # constant
    g = 9.8
    L = 0.1
    m = 1

    parser = argparse.ArgumentParser()
    parser.add_argument("-s",
                        "--step",
                        type=float,
                        help='the size time step, unit: second',
                        default=0.01)
    parser.add_argument("-t",
                        "--timelen",
                        type=int,
                        help='amount of time steps, unit: none',
                        default=100)
    parser.add_argument("-a",
                        "--angle",
                        type=int,
                        help='initial theta value, unit: degree',
                        default=6)
    args = parser.parse_args()

    # initial condition
    h = args.step
    time_length = args.timelen
    theta_0 = args.angle
    theta_0 = float(theta_0 * math.pi / 180)
    # exact energy
    E_exact = m * g * L * (1 - math.cos(theta_0))
    time, theta, time_half, E_half = compute_theta_energy(theta_0, h, time_length)
    make_plots(time, theta, time_half, E_half, E_exact, h, args.angle)
