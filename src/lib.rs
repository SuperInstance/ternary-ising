#![forbid(unsafe_code)]
//! Ternary Ising model — spins on a lattice where 0 is the topological insulator.

use std::collections::HashMap;

/// Ternary spin value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Spin { Down = -1, Zero = 0, Up = 1 }

impl Spin {
    pub fn val(self) -> i8 { self as i8 }
    pub fn from_i8(v: i8) -> Self { match v { -1 => Spin::Down, 0 => Spin::Zero, _ => Spin::Up } }
    pub fn is_insulator(self) -> bool { self == Spin::Zero }
    pub fn random(rng: &mut impl FnMut() -> f64) -> Self {
        let r = rng();
        if r < 0.333 { Spin::Down } else if r < 0.667 { Spin::Zero } else { Spin::Up }
    }
}

/// 2D ternary Ising lattice.
pub struct TernaryIsing {
    pub spins: Vec<Spin>,
    pub width: usize,
    pub height: usize,
    pub beta: f64,  // Inverse temperature
    pub j: f64,     // Coupling constant
}

impl TernaryIsing {
    pub fn new(width: usize, height: usize, beta: f64, j: f64) -> Self {
        let mut rng_state: u64 = 42;
        let mut rng = || -> f64 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };
        let spins = (0..width * height).map(|_| Spin::random(&mut rng)).collect();
        Self { spins, width, height, beta, j }
    }

    pub fn new_uniform(width: usize, height: usize, spin: Spin, beta: f64, j: f64) -> Self {
        let spins = vec![spin; width * height];
        Self { spins, width, height, beta, j }
    }

    fn idx(&self, r: usize, c: usize) -> usize { r * self.width + c }

    /// Get neighbors (4-connected, wrapping).
    pub fn neighbors(&self, r: usize, c: usize) -> Vec<(usize, usize)> {
        let h = self.height;
        let w = self.width;
        vec![
            ((r + h - 1) % h, c),
            ((r + 1) % h, c),
            (r, (c + w - 1) % w),
            (r, (c + 1) % w),
        ]
    }

    /// Local energy at site (r,c).
    pub fn local_energy(&self, r: usize, c: usize) -> f64 {
        let s = self.spins[self.idx(r, c)].val() as f64;
        let mut sum = 0.0;
        for (nr, nc) in self.neighbors(r, c) {
            let sn = self.spins[self.idx(nr, nc)].val() as f64;
            sum += sn;
        }
        -self.j * s * sum
    }

    /// Total energy of the system.
    pub fn total_energy(&self) -> f64 {
        let mut e = 0.0;
        for r in 0..self.height {
            for c in 0..self.width {
                let s = self.spins[self.idx(r, c)].val() as f64;
                // Only right and down to avoid double counting
                let right = self.spins[self.idx(r, (c + 1) % self.width)].val() as f64;
                let down = self.spins[self.idx((r + 1) % self.height, c)].val() as f64;
                e -= self.j * s * (right + down);
            }
        }
        e
    }

    /// Magnetization (signed).
    pub fn magnetization(&self) -> f64 {
        self.spins.iter().map(|s| s.val() as f64).sum::<f64>() / self.spins.len() as f64
    }

    /// Absolute magnetization.
    pub fn abs_magnetization(&self) -> f64 {
        self.spins.iter().map(|s| (s.val() as f64).abs()).sum::<f64>() / self.spins.len() as f64
    }

    /// Fraction of insulator (0) spins.
    pub fn insulator_fraction(&self) -> f64 {
        self.spins.iter().filter(|s| s.is_insulator()).count() as f64 / self.spins.len() as f64
    }

    /// Shannon entropy over {-1, 0, +1} distribution.
    pub fn entropy(&self) -> f64 {
        let mut counts: HashMap<Spin, usize> = HashMap::new();
        for &s in &self.spins { *counts.entry(s).or_insert(0) += 1; }
        let n = self.spins.len() as f64;
        let mut h = 0.0;
        for &c in counts.values() {
            let p = c as f64 / n;
            if p > 0.0 { h -= p * p.ln(); }
        }
        h
    }

    /// One Metropolis-Hastings sweep (N random flip attempts).
    pub fn metropolis_sweep(&mut self) -> usize {
        let n = self.spins.len();
        let mut accepted = 0;
        let mut rng_state: u64 = 7919;
        let mut rng = || -> f64 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        for _ in 0..n {
            let site = (rng() * n as f64) as usize;
            let r = site / self.width;
            let c = site % self.width;

            let old_spin = self.spins[site];
            let old_energy = self.local_energy(r, c);

            // Propose new spin (different from current)
            let new_spin = loop {
                let s = Spin::random(&mut rng);
                if s != old_spin { break s; }
            };

            self.spins[site] = new_spin;
            let new_energy = self.local_energy(r, c);

            let delta_e = new_energy - old_energy;
            if delta_e <= 0.0 || rng() < (-self.beta * delta_e).exp() {
                accepted += 1; // Accept
            } else {
                self.spins[site] = old_spin; // Reject
            }
        }
        accepted
    }

    /// Run N sweeps and return magnetization history.
    pub fn run(&mut self, sweeps: usize) -> Vec<f64> {
        let mut mag_history = Vec::with_capacity(sweeps);
        for _ in 0..sweeps {
            self.metropolis_sweep();
            mag_history.push(self.abs_magnetization());
        }
        mag_history
    }

    /// Finite-size scaling: estimate critical temperature by finding peak in susceptibility.
    pub fn estimate_tc(sizes: &[usize], beta_range: &[f64], sweeps_per_point: usize) -> Vec<(f64, f64)> {
        let mut susceptibility = Vec::new();
        for &beta in beta_range {
            let mut chi_sum = 0.0;
            for &size in sizes {
                let mut model = Self::new(size, size, beta, 1.0);
                let mags: Vec<f64> = model.run(sweeps_per_point);
                let mean_m = mags.iter().sum::<f64>() / mags.len() as f64;
                let var_m = mags.iter().map(|m| (m - mean_m).powi(2)).sum::<f64>() / mags.len() as f64;
                chi_sum += var_m * (size * size) as f64;
            }
            susceptibility.push((beta, chi_sum / sizes.len() as f64));
        }
        susceptibility
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spin_values() {
        assert_eq!(Spin::Down.val(), -1);
        assert_eq!(Spin::Zero.val(), 0);
        assert_eq!(Spin::Up.val(), 1);
        assert!(Spin::Zero.is_insulator());
    }

    #[test]
    fn test_uniform_lattice() {
        let model = TernaryIsing::new_uniform(5, 5, Spin::Up, 1.0, 1.0);
        assert!(model.spins.iter().all(|&s| s == Spin::Up));
        assert!((model.magnetization() - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_energy_ferromagnetic() {
        let model = TernaryIsing::new_uniform(4, 4, Spin::Up, 1.0, 1.0);
        let e = model.total_energy();
        // Each bond: -J * 1 * 1 = -1. Total bonds = 2*N = 32 (with wrapping)
        assert!(e < 0.0, "Ferromagnetic should have negative energy, got {}", e);
    }

    #[test]
    fn test_neighbors_wrapping() {
        let model = TernaryIsing::new(3, 3, 1.0, 1.0);
        let n = model.neighbors(0, 0);
        assert_eq!(n.len(), 4);
    }

    #[test]
    fn test_insulator_fraction() {
        let model = TernaryIsing::new_uniform(3, 3, Spin::Zero, 1.0, 1.0);
        assert!((model.insulator_fraction() - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_entropy_pure() {
        let model = TernaryIsing::new_uniform(4, 4, Spin::Up, 1.0, 1.0);
        assert!(model.entropy().abs() < 0.01, "Pure state should have zero entropy");
    }

    #[test]
    fn test_metropolis_sweep() {
        let mut model = TernaryIsing::new_uniform(5, 5, Spin::Up, 0.1, 1.0); // High temp
        let accepted = model.metropolis_sweep();
        assert!(accepted > 0, "Should accept some flips at high temperature");
    }

    #[test]
    fn test_low_temp_ordering() {
        let mut model = TernaryIsing::new_uniform(8, 8, Spin::Up, 5.0, 1.0); // Low temp
        model.run(50);
        // At low temperature, should stay mostly ordered
        assert!(model.abs_magnetization() > 0.5, "Should stay ordered at low T");
    }

    #[test]
    fn test_high_temp_disorder() {
        let mut model = TernaryIsing::new(8, 8, 0.5, 1.0); // High temp
        model.run(200);
        let m = model.abs_magnetization();
        assert!(m < 0.9, "Should not be fully ordered at high T, got m={}", m);
    }

    #[test]
    fn test_local_energy() {
        let model = TernaryIsing::new_uniform(3, 3, Spin::Up, 1.0, 1.0);
        let e = model.local_energy(1, 1);
        assert!((e - (-4.0)).abs() < 0.01, "Center of uniform Up should be -4J, got {}", e);
    }

    #[test]
    fn test_entropy_mixed() {
        let mut model = TernaryIsing::new(10, 10, 100.0, 1.0);
        model.run(200);
        let h = model.entropy();
        assert!(h > 0.5, "Mixed state should have significant entropy, got {}", h);
    }

    #[test]
    fn test_run_returns_history() {
        let mut model = TernaryIsing::new_uniform(4, 4, Spin::Up, 1.0, 1.0);
        let history = model.run(10);
        assert_eq!(history.len(), 10);
        assert!(history.iter().all(|&m| m >= 0.0 && m <= 1.0));
    }

    #[test]
    fn test_zero_spin_contribution() {
        let mut model = TernaryIsing::new_uniform(3, 3, Spin::Zero, 1.0, 1.0);
        let e = model.total_energy();
        assert!(e.abs() < 0.01, "All-zero should have zero energy, got {}", e);
    }
}
