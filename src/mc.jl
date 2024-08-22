using LoopVectorization

function metropolis(nx, ny, beta, spins)
    for i = 1:ny
      for j = (1 + (i & 1)):2:nx
        local_update(nx, ny, j, i, beta, spins)
      end
    end
    for i = 1:ny
      for j = (2 - (i & 1)):2:nx
        local_update(nx, ny, j, i, beta, spins)
      end
    end
end

function local_update(nx, ny, x, y, beta, spins)
    dy = [1, 0, -1, 0]
    dx = [0, 1, 0, -1]
    near_summ = 0
    for d = 1:4
        near_y = y + dy[d]
        if near_y > ny
            near_y = 1
        elseif near_y < 1
            near_y = ny
        end
        near_x = x + dx[d]
        if near_x > nx
            near_x = 1
        elseif near_x < 1
            near_x = nx
        end
        near_summ += spins[near_x, near_y]
    end
    energy_diff = 2 * spins[x, y] * near_summ
    if (energy_diff <= 0)
        spins[x, y] = - spins[x, y]
    elseif (rand() < exp(-beta * energy_diff))
        spins[x, y] = - spins[x, y]
    end
end

function calc_magne(nx, ny, nall_inv, spins)
    summ = 0
    for y = 1:ny
        for x = 1:nx
            summ += spins[x][y]
        end
    end
    summ * nall_inv
end

function main(nx, ny, nall_inv, beta, mcs, nsample)
    magnes = zeros(Float64, mcs)
    for j = 1:nsample
        spins = ones(Int32, nx, ny)

        for i = 1:mcs
            metropolis(nx, ny, beta, spins)
            magnes[i] += calc_magne(nx, ny, nall_inv, spins)
        end
    end
    for i = 1:mcs
        println(i, magnes[i])
    end
    # println(spins)
end

const nx = 100
const ny = ny
const nall = nx * ny
const nall_inv = 1 / nall
const kbt = 2.3
const beta = 1 / kbt
const mcs = 1000
const nsample = 100
@time main(nx, ny, nall_inv, beta, mcs, nsample)
