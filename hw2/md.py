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
import matplotlib.pyplot as plt

rho = float(1458 / 1680)  # density of Argon in reduced units
dt = 0.004  # time step size
T_0 = float(203 / 119.8)  # temperature
N_cell = 3  # number of fcc unitcells in one direction
N = 4 * N_cell**3  # the total number of particles in the system
L_box = (N / rho)**(1 / 3.0)  # length of the whole simulation box
L_cell = L_box / N_cell  # length of a unitcell
F = np.zeros((N, N, 3))  # matrix that contains all forces
ind = np.triu_indices(N, k=1)  # indices of upper triangular matrix
r_final = np.zeros((N, N))


def IC_pos(N_cell, L_cell):
    '''
    use fcc structure to initilize positions
    '''
    pos = [[[x, y, z], [x, 0.5 + y, 0.5 + z], [0.5 + x, y, 0.5 + z], [0.5 + x, 0.5 + y, z]]
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
    # print("before asign the F shape")
    # print(F.shape)
    r_vec = pos[ind[0]] - pos[ind[1]]
    # why here?
    r_vec = r_vec - np.rint(r_vec / L_box) * L_box
    r_sq = np.sum(r_vec**2, axis=1)
    F_vec = -(48 / r_sq**7 - 24 / r_sq**4)[:, None] * r_vec
    r_final[ind[0], ind[1]] = r_sq**(1 / 2)

    F[ind[0], ind[1]] = F_vec
    # print("in find force, the F shape")
    # print(F.shape)
    # print(id(F))
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

    # print("in time step, the F shape is")
    # print(F.shape)
    # print(id(F))
    vel += 0.5 * F * dt
    kin = 0.5 * np.sum(vel**2)
    # potential and kenetic energy
    return pos, vel, F, pot, kin, P, pos_folded


def simulate(file):
    kins, pots, Ps, msd_list = [], [], [], []
    vac = 0
    pos = IC_pos(N_cell, L_cell)
    pos_0 = pos
    vel = IC_vel(N)
    vel_0 = vel
    F = find_force(pos)[0]
    for t in range(2000):
        pos, vel, F, pot, kin, P, pos_folded = time_step(pos, vel, F)
        #print(type(vel_0), vel_0.shape)
        #print(type(vel), vel.shape)
        vac += np.sum(vel_0 * vel) * dt
        if t > 1000:  # production run
            #print("positon shape", pos.shape)
            #print("vel shape", vel.shape)
            file.write("ARGON MD, t=   %.5f\n" % (t * dt))
            file.write("  108\n")
            for i in range(N):
                file.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" %
                           (i, 'ARGON', 'A', i, pos_folded[i][0], pos_folded[i][1],
                            pos_folded[i][2], vel[i][0], vel[i][1], vel[i][2]))
            file.write("%8.5f%8.5f%8.5f\n" % (L_box, L_box, L_box))
            msd = np.sum((pos - pos_0)**2) / N
            msd_list.append(msd)
            kins.append(kin)
            pots.append(pot)
            Ps.append(P)
        else:  # equillirum run
            vel *= np.sqrt(N * 3 * T_0 / (2 * kin))
    return np.array(kins), np.array(pots), np.array(Ps), np.array(msd_list), vac


def radial_distribution(r_final):
    # calculate the radial distribution of the first atom
    r_distribution = (r_final + np.transpose(r_final)).reshape(-1)
    r_ = np.linspace(np.min(r_distribution), np.max(r_distribution), 100)
    dr = (np.max(r_distribution) - np.min(r_distribution)) / 100
    v = r_**2 * (4 * math.pi * dr)
    cout = np.zeros(100)

    print(r_distribution)
    for r in r_distribution:
        for i in range(99):
            if (r_[i] < r < r_[i + 1]):
                cout[i] += 1

    return np.array(cout) / (N * np.array(v)), r_


# The simulation starts here
if __name__ == "__main__":
    fo = open("./argon.gro", "w")
    kins, pots, Ps, msds, vacs = simulate(fo)
    fo.close()
    totol_energy = pots + kins
    times = np.arange(1000, 1999) * 0.004
    if (times.shape[0] != totol_energy.shape[0]):
        print(totol_energy.shape)
        print(times.shape)
        print("wrong")
    else:
        plt.plot(times, pots, label="potential energy")
        plt.plot(times, kins, label="kinetic engergy")
        plt.plot(times, totol_energy, label="total energy")
        plt.xlabel("Time: reduced unit")
        plt.ylabel("Energy: reduced unit $\epsilon$")
        plt.legend(loc='upper right')
        plt.savefig("./figs/energy_time.png", dpi=400)
        plt.clf()

    couts, r_dis = radial_distribution(r_final)
    if (couts.shape[0] == r_dis.shape[0]):
        plt.plot(r_dis, couts)
        plt.xlabel("$\\rho(\sigma)$")
        plt.ylabel("$g(r)$")
        plt.savefig("./figs/raditial_func.png", dpi=400)
        plt.clf()

    D1 = np.average(np.gradient(msds)/6)/0.004
    D2 = vacs / (3 * N)
    print(D1, D2)

    T = np.mean(kins * 2 / (3 * N))  # temperature
    P = 1 - np.mean(Ps) / (3 * N * T) - 16 * np.pi * rho / (3 * T * L_box**3
                                                            )  # compressibility factor
    P = P * T * rho  # pressure here

    #print(T, P)  # how about thermal fluctuation?
    #print(pots)
