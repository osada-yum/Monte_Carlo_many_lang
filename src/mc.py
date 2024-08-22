import sys

import numpy as np

def metropolis(nx, ny, beta, spins):
    for i in range(nx):
        for j in range(0 + (i & 1), ny, 2):
            local_flip(nx, ny, i, j, beta, spins)
    for i in range(nx):
        for j in range(1 - (i & 1), ny, 2):
            local_flip(nx, ny, i, j, beta, spins)

def local_flip(nx, ny, x, y, beta, spins):
    dx = np.array([1, 0, -1, 0])
    dy = np.array([0, 1, 0, -1])
    near_summ = 0
    for d in range(4):
        near_x = x + dx[d]
        if near_x >= nx:
            near_x = 0
        elif near_x < 0:
            near_x = nx - 1
        near_y = y + dy[d]
        if near_y >= ny:
            near_y = 0
        elif near_y < 0:
            near_y = ny - 1
        near_summ += spins[near_x][near_y]
    delta_energy = 2 * spins[x][y] * near_summ
    if delta_energy <= 0:
        spins[x][y] = - spins[x][y]
    elif np.exp(- beta * delta_energy) <= np.random.rand():
        spins[x][y] = - spins[x][y]

def calc_magne(nx, ny, nall_inv, spins):
    return np.sum(spins, axis = (0, 1)) * nall_inv


nx = 100
ny = nx
nall = nx * ny
nall_inv = 1 / nall

kbt = 2.3
beta = 1 / kbt

spins = np.ones(nx * ny, dtype = np.int32).reshape(nx, ny)

mcs = 1000
nsample = 100

magnes = np.zeros(mcs, dtype = np.float64)

for j in range(nsample):
    print(f"sample {j + 1}", file = sys.stderr)
    for i in range(mcs):
        metropolis(nx, ny, beta, spins)
        magnes[i] += calc_magne(nx, ny, nall_inv, spins)

for i in range(mcs):
    print(f"{i + 1} {magnes[i]}")
# print(spins)
