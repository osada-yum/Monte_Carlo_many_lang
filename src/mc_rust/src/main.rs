use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

const NX: usize = 100_usize;
const NY: usize = 100_usize;
const NALL: usize = NX * NY;
const NALL_INV: f64 = 1_f64 / (NALL as f64);

const KBT: f64 = 2.3_f64;
const BETA: f64 = 1_f64 / KBT;

const ND: usize = 4_usize;
const DX: [isize; ND] = [1, 0, -1, 0];
const DY: [isize; ND] = [0, 1, 0, -1];

fn metropolis(spins: &mut Vec<Vec<i32>>, rng: &mut ChaCha8Rng) {
    for i in 0..NX {
        let start = i & 1_usize;
        for j in (start..NY).step_by(2) {
            local_flip(i, j, spins, rng);
        }
    }
    for i in 0..NX {
        let start = (i + 1) & 1_usize;
        for j in (start..NY).step_by(2) {
            local_flip(i, j, spins, rng);
        }
    }
}

fn local_flip(x: usize, y: usize, spins: &mut Vec<Vec<i32>>, rng: &mut ChaCha8Rng) {
    let mut near_summ = 0;
    let cx = x as isize;
    let cy = y as isize;
    for d in 0..ND {
        let nx = cx + DX[d];
        let ny = cy + DY[d];
        let nx: usize = if nx < 0 {
            NX - 1
        } else if nx >= NX as isize {
            0
        } else {
            nx as usize
        };
        let ny: usize = if ny < 0 {
            NY - 1
        } else if ny >= NY as isize {
            0
        } else {
            ny as usize
        };
        near_summ += spins[nx][ny];
    }
    let delta_energy = 2 * spins[x][y] * near_summ;
    if delta_energy <= 0 {
        spins[x][y] = -spins[x][y]
    } else {
        let r: f64 = rng.random();
        if r < (-BETA * (delta_energy as f64)).exp() {
            spins[x][y] = -spins[x][y];
        }
    }
}

fn calc_magne(spins: &Vec<Vec<i32>>) -> f64 {
    let mut summ = 0;
    for i in 0..NX {
        for j in 0..NY {
            summ += spins[i][j];
        }
    }
    (summ as f64) * NALL_INV
}

fn main() {
    const MCS: usize = 1000;
    const NSAMPLE: usize = 100;

    let mut spins: Vec<Vec<i32>>;
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(42);

    let mut magnes: Vec<f64> = vec![0_f64; MCS];

    for s in 0..NSAMPLE {
        eprintln!("sample: {}", s + 1);
        spins = vec![vec![1; NY]; NX];
        for t in 0..MCS {
            // eprintln!("MCS {}...", t + 1);
            metropolis(&mut spins, &mut rng);
            magnes[t] += calc_magne(&spins);
        }
    }
    for t in 0..MCS {
        magnes[t] = magnes[t] / (NSAMPLE as f64);
        println!("{} {}", t, magnes[t]);
    }
}
