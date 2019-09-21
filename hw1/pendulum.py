import numpy as np
import matplotlib.pyplot as plt
import math
"""
The equation we want to solve:
\theta(t+h)=2\theta(t)-\theta(t-h)+(-g/l)\theta(t)h^2
\theta(t=0)=1/30
"""

# constant
g = 9.8 
L = 0.1
h = 0.01
m = 1

# initial condition
theta_0 = float(1 / 30)


def theta_generator(steps, small_angle=True):
    # initial condition
    i, theta_t_minus_h, theta_t = 0, theta_0, theta_0
    while i < steps:
        yield theta_t
        # Verlet method
        if small_angle:
            theta_t_plus_h = 2 * theta_t - theta_t_minus_h - g / L * theta_t * h**2
        else:
            theta_t_plus_h = 2 * theta_t - theta_t_minus_h - g / L * math.sin(theta_t) * h**2
        # propagte
        i, theta_t, theta_t_minus_h = i + 1, theta_t_plus_h, theta_t


# theta list with time
theta = [x for x in theta_generator(100)]
time = [t for t in range(100)]

# theta + h/2 with time, half represent time at t + h/2
time_half = [t + h / 2 for t in time]
time_half.pop(-1)
theta_half = [(theta[i + 1] - theta[i]) / 2 for i in range(len(theta))]
velocity_half = [L * (theta[i + 1] - theta[i]) / (2 * h) for i in range(len(theta))]
T_half = [1 / 2 * m * v**2 for v in velocity_half]
V_half = [m * g * L * (1 - math.cos(theta)) for theta in theta_half]
E_half = np.array(T_half) + np.array(V_half)
assert len(T_half)==len(V_half), "Error"

print(len(theta))
print(len(time))
plt.plot(time, theta)