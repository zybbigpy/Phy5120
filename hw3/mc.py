'''
Monte Carlo Simulation of london jones fluid
Reduced unit: \sigma=\epsilon=mass=k_B=1
'''

import math
import numpy as np
from itertools import product
import matplotlib.pyplot as plt

rho = float(1458 / 1680)  # density of Argon in reduced units
T_0 = float(203 / 119.8)  # temprature
N_cell = 6  # number of fcc unitcells in one direction
N = 4 * N_cell**3  # the total number of particles in the system
L_box = (N / rho)**(1 / 3.0)  # length of the whole simulation box
L_cell = L_box / N_cell  # length of a unitcell
ind = np.triu_indices(N, k=1)  # indices of upper triangular matrix
beta = 1 / T_0


def IC_pos(N_cell, L_cell):
    '''
    use fcc structure to initilize positions
    '''
    pos = [[[x, y, z], [x, 0.5 + y, 0.5 + z], [0.5 + x, y, 0.5 + z], [0.5 + x, 0.5 + y, z]]
           for x, y, z in product(range(N_cell), range(N_cell), range(N_cell))]
    pos = np.array(pos).reshape((-1, 3))

    return pos * L_cell


def find_potential(pos, L_box=L_box):
    '''
    find the whole potential of the simulating box
    '''

    r_vec = pos[ind[0]] - pos[ind[1]]
    r_vec = r_vec - np.rint(r_vec / L_box) * L_box
    r_sq = np.sum(r_vec**2, axis=1)

    pot = np.sum(4 / r_sq**6 - 4 / r_sq**3)

    return pot


def time_step(pos):
    '''
    use monte carlo method instead of md
    '''

    pos_folded = np.mod(pos, L_box)
    E_o = find_potential(pos_folded)

    pos_change = (np.random.rand(pos.shape[0], pos.shape[1]) - 0.5) * 0.001
    pos_n = pos + pos_change

    posn_folded = np.mod(pos_n, L_box)
    E_n = find_potential(posn_folded)

    # accept possibility
    # print(-beta*(E_n-E_o))
    prob = math.exp(-beta * (E_n - E_o))

    acc = min(1, prob)
    randnum = np.random.rand()
    #print("E_n[%.3f] E_o[%.3f] prob[%.3f] randnum[%.3f], acc[%.3f]"%(E_n, E_o, prob, randnum, acc))

    if (acc >= randnum):
        pos += pos_change

    return pos


def get_rfinal(pos):
    '''
    calculate r final
    '''

    r_vec = pos[ind[0]] - pos[ind[1]]
    r_vec = r_vec - np.rint(r_vec / L_box) * L_box

    r_sq = np.sum(r_vec**2, axis=1)
    r_final = np.zeros((N, N))
    r_final[ind[0], ind[1]] = r_sq**(1 / 2)

    return r_final


def simulate(steps):
    '''
    Monte Calor Simulation
    '''
    pos = IC_pos(N_cell, L_cell)
    for t in range(steps):
        print(t)
        pos = time_step(pos)

    return get_rfinal(pos)


def radial_distribution(r_final):
    '''
    calculate the radial distribution of the first atom
    '''

    r_distribution = (r_final + np.transpose(r_final)).reshape(-1)
    r_ = np.linspace(np.min(r_distribution), np.max(r_distribution), 100)
    dr = (np.max(r_distribution) - np.min(r_distribution)) / 100
    v = (r_**2 * (4 * math.pi * dr))[1:]
    cout = np.zeros(99)

    for r in r_distribution:
        for i in range(98):
            if (r_[i] < r < r_[i + 1]):
                cout[i] += 1

    return np.array(cout) / (N * np.array(v)), r_[1:]


if __name__ == "__main__":

    steps = 2000
    r_final = simulate(steps)

    couts, r_dis = radial_distribution(r_final)

    # print("L_cell[%.3f], L_box[%.3f]"%(L_cell, L_box))
    print(couts.shape[0], r_dis.shape[0])
    if (couts.shape[0] == r_dis.shape[0]):
        print("i am here")
        plt.plot(r_dis[:65], couts[:65])
        plt.xlabel("$\\rho(\sigma)$")
        plt.ylabel("$g(r)$")
        plt.savefig("./figs/raditial_func.png", dpi=400)
        plt.clf()
