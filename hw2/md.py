'''
MD simulation of Lennard-Jones fluids
modified from
https://github.com/basnijholt/molecular-dynamics-Python/blob/master/MD.ipynb
The original code is wrong. If you simply copy it, you will get ZERO mark!
Reduced unit: \sigma=\epsilon=mass=k_B=1
'''

import numpy as np
from itertools import product
rho = 0.6  # density of Argon in reduced units
dt = 0.004  # time step size
T_0 = 2  # temperature
N_cell = 3  # number of fcc unitcells in one direction
N = 4 * N_cell**3  # the total number of particles in the system
L_box = (N / rho)**(1 / 3.0)  # length of the whole simulation box
L_cell = L_box / N_cell  # length of a unitcell
F = np.zeros((N, N, 3))  # matrix that contains all forces
ind = np.triu_indices(N, k=1)  # indices of upper triangular matrix


def IC_pos(N_cell, L_cell):
    '''
    use fcc structure to initilize positions
    '''
    pos = [[[x, y, z], [x, 0.5 + y, 0.5 + z], [0.5 + x, y, 0.5 + z],
            [0.5 + x, 0.5 + y, z]]
           for x, y, z in product(range(N_cell), range(N_cell), range(N_cell))]
    pos = np.array(pos).reshape((-1, 3))
    return pos * L_cell


def IC_vel(N):
    '''
    Maxwell-Boltzman distribution is a normal distribution
    '''
    vel = np.sqrt(T_0) * np.random.randn(N, 3)
    vel -= np.average(vel, axis=0)
    return vel


def find_force(pos, L_box=L_box):
    '''
    Minimum image convention. 
    '''
    r_vec = pos[ind[0]] - pos[ind[1]]
    r_vec = r_vec - np.rint(r_vec / L_box) * L_box
    r_sq = np.sum(r_vec**2, axis=1)
    F_vec = -(48 / r_sq**7 - 24 / r_sq**4)[:, None] * r_vec
    F[ind[0], ind[1]] = F_vec
    pot = np.sum(4 / r_sq**6 - 4 / r_sq**3)
    P = np.sum(F_vec * r_vec)
    return np.sum(F, axis=0) - np.sum(F, axis=1), pot, P


def time_step(pos, vel, F):
    vel += 0.5 * F * dt
    pos = pos + vel * dt
    pos_folded = np.mod(pos, L_box)
    #pos = np.mod(pos + vel * dt, L_box) # why both pos and pos_folded?
    F, pot, P = find_force(pos_folded)

    vel += 0.5 * F * dt
    kin = 0.5 * np.sum(vel**2)
    return pos, vel, F, pot, kin, P


def simulate():
    kins, pots, Ps = [], [], []
    pos = IC_pos(N_cell, L_cell)
    vel = IC_vel(N)
    F = find_force(pos)[0]
    for t in range(2000):
        pos, vel, F, pot, kin, P = time_step(pos, vel, F)
        if t > 1000:  # production run
            kins.append(kin)
            pots.append(pot)
            Ps.append(P)
        else:  # equillirum run
            vel *= np.sqrt(N * 3 * T_0 / (2 * kin))
    return np.array(kins), np.array(pots), np.array(Ps)


# The simulation starts here
if __name__ == "__main__":
    kins, pots, Ps = simulate()

    T = np.mean(kins * 2 / (3 * N))  # temperature
    P = 1 - np.mean(Ps) / (3 * N * T) - 16 * np.pi * rho / (
        3 * T * L_box**3)  # compressibility factor
    P = P * T * rho  # pressure here

    print(T, P)  # how about thermal fluctuation?
    print(pots)
