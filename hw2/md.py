'''
MD simulation of Lennard-Jones fluids
modified from
https://github.com/basnijholt/molecular-dynamics-Python/blob/master/MD.ipynb
The original code is wrong. If you simply copy it, you will get ZERO mark!
Reduced unit: \sigma=\epsilon=mass=k_B=1
'''

import numpy as np
from itertools import product
import math


rho = 0.6  # density of Argon in reduced units
dt = 0.004  # time step size
T_0 = 2  # temperature
N_cell = 3  # number of fcc unitcells in one direction
N = 4 * N_cell**3  # the total number of particles in the system
L_box = (N / rho)**(1 / 3.0)  # length of the whole simulation box
L_cell = L_box / N_cell  # length of a unitcell
F = np.zeros((N, N, 3))  # matrix that contains all forces
ind = np.triu_indices(N, k=1)  # indices of upper triangular matrix
r_final = np.zeros((N,N))


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
    print("before asign the F shape")
    print(F.shape)
    r_vec = pos[ind[0]] - pos[ind[1]]
    # why here?
    r_vec = r_vec - np.rint(r_vec / L_box) * L_box
    r_sq = np.sum(r_vec**2, axis=1)
    F_vec = -(48 / r_sq**7 - 24 / r_sq**4)[:, None] * r_vec
    r_final[ind[0], ind[1]] = r_sq**(1/2)

    F[ind[0], ind[1]] = F_vec
    print("in find force, the F shape")
    print(F.shape)
    print(id(F))
    pot = np.sum(4 / r_sq**6 - 4 / r_sq**3)
    P = np.sum(F_vec * r_vec)
    #  why minus instead of plus
    return np.sum(F, axis=0) - np.sum(F, axis=1), pot, P


def time_step(pos, vel, F):
    vel += 0.5 * F * dt
    pos = pos + vel * dt
    pos_folded = np.mod(pos, L_box)
    #pos = np.mod(pos + vel * dt, L_box) # why both pos and pos_folded?
    F, pot, P = find_force(pos_folded)

    print("in time step, the F shape is")
    print(F.shape)
    print(id(F))
    vel += 0.5 * F * dt
    kin = 0.5 * np.sum(vel**2)
    # potential and kenetic energy
    return pos, vel, F, pot, kin, P


def simulate():
    kins, pots, Ps, msd_list = [], [], [], []
    pos = IC_pos(N_cell, L_cell)
    pos_0 = pos
    vel = IC_vel(N)
    F = find_force(pos)[0]
    for t in range(2000):
        pos, vel, F, pot, kin, P = time_step(pos, vel, F)
        if t > 1000:  # production run
            msd = np.sum((pos-pos_0)**2)/N
            msd_list.append(msd)
            kins.append(kin)
            pots.append(pot)
            Ps.append(P)
        else:  # equillirum run
            vel *= np.sqrt(N * 3 * T_0 / (2 * kin))
    return np.array(kins), np.array(pots), np.array(Ps), np.array(msd_list)

def radial_distribution(r_final):
    # calculate the radial distribution of the first atom
    r_distribution = np.array(r_final[0])
    r_ = np.linspace(1, np.max(r_distribution), 100)
    dr = (np.max(r_distribution) - 1)/100
    v = r_**2*4*math.pi*dr
    cout = np.zeros(100)

    for r in r_distribution:
        for i in range(100):
            if(r>r_[i]):
                cout[i] +=1

    return np.array(cout)/np.array(v), r_

# The simulation starts here
if __name__ == "__main__":
    kins, pots, Ps, msds = simulate()

    T = np.mean(kins * 2 / (3 * N))  # temperature
    P = 1 - np.mean(Ps) / (3 * N * T) - 16 * np.pi * rho / (
        3 * T * L_box**3)  # compressibility factor
    P = P * T * rho  # pressure here

    print(T, P)  # how about thermal fluctuation?
    print(pots)
