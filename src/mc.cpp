#include <iostream>
#include <iomanip>
#include <cstdint>
#include <vector>
#include <cmath>
#include <random>

constexpr const std::int64_t NX = 100LL, NY = NX, NALL = NX * NY;
constexpr const double NALL_INV = 1.0 / NALL;
constexpr const double KBT = 2.3, BETA = 1.0 / KBT;

constexpr const std::int32_t ND = 4;
constexpr const std::int32_t DX[ND] = {1, 0, -1, 0};
constexpr const std::int32_t DY[ND] = {0, 1, 0, -1};

using Lattice = std::vector<std::vector<std::int32_t> >;

void set_order_spin(Lattice& spins) {
  for (std::int64_t x = 0; x < NX; ++x) {
    for (std::int64_t y = 0; y < NX; ++y) {
      spins[x][y] = 1;
    }
  }
}

void local_update(Lattice& spins, const std::int64_t x, const std::int64_t y, std::mt19937& mt, std::uniform_real_distribution<double>& dist) {
  std::int32_t near_summ = 0;
  for (std::int32_t d = 0; d < ND; ++d) {
    std::int32_t nx = x + DX[d];
    if (nx < 0) {
      nx = NX - 1;
    } else if (nx >= NX) {
      nx = 0;
    }
    std::int32_t ny = y + DY[d];
    if (ny < 0) {
      ny = NY - 1;
    } else if (ny >= NY) {
      ny = 0;
    }
    near_summ += spins[nx][ny];
  }
  const std::int32_t delta_energy = 2 * spins[x][y] * near_summ;
  if (delta_energy <= 0) {
    spins[x][y] = - spins[x][y];
  } else {
    const double r = dist(mt);
    if (r < std::exp(- BETA * delta_energy)) {
      spins[x][y] = - spins[x][y];
    }
  }
}

void update(Lattice& spins, std::mt19937& mt, std::uniform_real_distribution<double>& dist) {
  for (std::int64_t i = 0LL; i < NX; ++i) {
    const std::int64_t start = i & 1LL;
    for (std::int64_t j = start; j < NY; j += 2LL) {
      local_update(spins, i, j, mt, dist);
    }
  }
  for (std::int64_t i = 0LL; i < NX; ++i) {
    const std::int64_t start = (i + 1LL) & 1LL;
    for (std::int64_t j = start; j < NY; j += 2LL) {
      local_update(spins, i, j, mt, dist);
    }
  }
}

double calc_magne(const Lattice& spins) {
  std::int64_t summ = 0LL;
  for (std::int64_t x = 0; x < NX; ++x) {
    for (std::int64_t y = 0; y < NX; ++y) {
      summ += spins[x][y];
    }
  }
  return summ * NALL_INV;
}

int main()
{
  constexpr const std::int32_t MCS = 1000;
  constexpr const std::int32_t NSAMPLE = 100;

  Lattice spins(NX, std::vector<std::int32_t>(NY));
  std::mt19937 mt(42);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<double> magnes(MCS, 0.0);

  for (std::int32_t i = 0; i < NSAMPLE; ++i) {
    set_order_spin(spins);
    std::cerr << "sample " << i + 1 << "\n";
    for (std::int32_t j = 0; j < MCS; ++j) {
      update(spins, mt, dist);
      magnes[j] += calc_magne(spins);
    }
  }
  for (std::int32_t j = 0; j < MCS; ++j) {
    magnes[j] /= NSAMPLE;
    std::cout << j << " " << std::setprecision(15) << magnes[j] << "\n";
  }
  return 0;
}
