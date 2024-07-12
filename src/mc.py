import numpy as np

def metropolis(nx, ny, beta, spins):
    for i in range(ny):
        for j in range(0 + (i & 1), nx, 2):
            local_flip(nx, ny, j, i, beta, spins)
    for i in range(ny):
        for j in range(1 - (i & 1), nx, 2):
            local_flip(nx, ny, j, i, beta, spins)

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


nx = 100
ny = 100
kbt = 2.2
beta = 1 / kbt
spins = np.ones(nx * ny).reshape(nx, ny)
mcs = 100

for i in range(mcs):
    metropolis(nx, ny, beta, spins)

print(spins)
