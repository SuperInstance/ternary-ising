# ternary-ising

Ising model simulation for **ternary spin systems** {-1, 0, +1} with Metropolis-Hastings dynamics and critical temperature detection. The zero spin state acts as a **topological insulator** — a non-magnetic phase that breaks the binary up/down symmetry of classical Ising models.

## Why It Matters

The classical Ising model (Ising, 1925; Onsager, 1944) has two spin states: up (+1) and down (−1). Adding a third state (0) creates a **Blume-Capel model** (Blume, 1966; Capel, 1966) with rich phase behavior:

- **Ferromagnetic phase**: Spins align (mostly +1 or mostly −1)
- **Paramagnetic phase**: Random spins (uniform mix of −1, 0, +1)
- **Insulator phase**: Dominated by 0-spins — no magnetic ordering
- **Triple critical point**: Where all three phases meet

The 0-spin is not just "no measurement" — it's a **topological insulator state** that interrupts magnetic interactions. Two neighboring +1 spins separated by a 0-spin do not interact, creating a lattice with position-dependent connectivity.

Applications include:
- Phase transitions in ternary magnetic materials
- Ternary neural network energy landscapes (weights {-1, 0, +1})
- Quantum annealing with three-state systems
- Social dynamics with "neutral" agents

## How It Works

### Hamiltonian

The ternary Ising Hamiltonian on a 2D lattice with periodic boundaries:

```
H = −J · Σ_{<i,j>} sᵢ · sⱼ
```

where sᵢ ∈ {−1, 0, +1}, J is the coupling constant, and ⟨i,j⟩ denotes nearest neighbors (4-connected).

**Local energy** at site (r,c):

```
E(r,c) = −J · s_{r,c} · Σ_{neighbors} s_{nᵣ,n_c}
```

**Total energy** (avoiding double-counting):

```
E = −J · Σᵢ sᵢ · (s_{right(i)} + s_{down(i)})
```

### Key Observables

**Magnetization** (signed):

```
m = (1/N) · Σᵢ sᵢ      ∈ [−1, 1]
```

**Absolute magnetization** (order parameter):

```
|m| = (1/N) · Σᵢ |sᵢ|   ∈ [0, 1]
```

Note: |m| measures magnetic order, not insulator fraction. A lattice of all 0-spins has |m| = 0.

**Insulator fraction**:

```
f₀ = (1/N) · Σᵢ 𝟙[sᵢ = 0]
```

**Shannon entropy** over the spin distribution:

```
S = −Σₛ p(s) · ln p(s)      where s ∈ {−1, 0, +1}
```

Maximum entropy S_max = ln(3) ≈ 1.0986 when all three states are equally populated.

### Metropolis-Hastings Dynamics

Each sweep attempts N = width × height random spin flips:

```
1. Pick random site (r, c)
2. Pick random new spin s' ∈ {−1, 0, +1} \ {s_current}
3. Compute ΔE = E_new − E_old = −J · (s' − s) · Σ_{neighbors} sₙ
4. Accept if ΔE < 0  OR  ξ < exp(−β · ΔE)   where ξ ~ Uniform(0,1)
5. Otherwise reject
```

where β = 1/(k_B T) is the inverse temperature.

**Acceptance probability**:

```
P(accept) = min(1, exp(−β · ΔE))
```

This satisfies **detailed balance**: P(s → s') · P(s) = P(s' → s) · P(s'), ensuring convergence to the Boltzmann distribution.

### Critical Temperature

For the classical 2D Ising model (binary spins), Onsager's exact solution gives:

```
T_c = 2J / (k_B · ln(1 + √2)) ≈ 2.269 J/k_B
```

For the ternary (Blume-Capel) model, T_c depends on the crystal field (energy cost of the 0-state). The crate's `metropolis_sweep` allows measuring observables across temperatures to locate T_c numerically.

### Complexity

| Operation | Time | Space |
|-----------|------|-------|
| `new(w, h, β, J)` | O(w·h) | O(w·h) |
| `local_energy(r, c)` | O(1) | O(1) |
| `total_energy()` | O(w·h) | O(1) |
| `magnetization()` | O(w·h) | O(1) |
| `entropy()` | O(w·h) | O(1) |
| `metropolis_sweep()` | O(w·h) | O(1) |

Each sweep is O(N) — one random flip attempt per spin, O(1) energy delta computation.

## Quick Start

```rust
use ternary_ising::{TernaryIsing, Spin};

// Create a 32×32 lattice at β = 0.5 (above T_c), J = 1.0
let mut lattice = TernaryIsing::new(32, 32, 0.5, 1.0);

println!("Initial energy: {:.2}", lattice.total_energy());
println!("Magnetization: {:.3}", lattice.magnetization());
println!("Insulator fraction: {:.3}", lattice.insulator_fraction());

// Run 1000 Metropolis sweeps to reach equilibrium
for _ in 0..1000 {
    let accepted = lattice.metropolis_sweep();
}

println!("Equilibrium energy: {:.2}", lattice.total_energy());
println!("Equilibrium magnetization: {:.3}", lattice.magnetization());
println!("Entropy: {:.3}", lattice.entropy());
```

### Temperature Sweep (Finding T_c)

```rust
// Scan temperature to find the phase transition
for beta in (1..100).map(|t| t as f64 / 50.0) {
    let mut lat = TernaryIsing::new(32, 32, beta, 1.0);
    for _ in 0..500 { lat.metropolis_sweep(); } // equilibrate
    println!("β={:.2}, |m|={:.3}, f₀={:.3}", beta, lat.abs_magnetization(), lat.insulator_fraction());
}
```

## API

### `TernaryIsing`

| Method | Description |
|--------|-------------|
| `new(w, h, beta, j)` | Random initial configuration |
| `new_uniform(w, h, spin, beta, j)` | Uniform initial state |
| `local_energy(r, c) -> f64` | Energy at site |
| `total_energy() -> f64` | Total system energy |
| `magnetization() -> f64` | Signed average spin |
| `abs_magnetization() -> f64` | Absolute average spin |
| `insulator_fraction() -> f64` | Fraction of 0-spins |
| `entropy() -> f64` | Shannon entropy |
| `metropolis_sweep() -> usize` | N random flip attempts, returns accept count |

### `Spin`

| Variant | Value | Semantic |
|---------|-------|----------|
| `Down` | −1 | Anti-aligned |
| `Zero` | 0 | Topological insulator |
| `Up` | +1 | Aligned |

## Architecture Notes

This crate implements **η (eta) layer** simulation in the γ + η = C framework:

- **η (eta)**: The Monte Carlo simulation engine — energy computation, spin dynamics, observable measurement. This crate provides the η-layer physics.
- **γ (gamma)**: External coordination — parallel tempering, replica exchange, and multi-lattice synchronization would be provided by ecosystem coordination crates.
- **C**: The complete statistical mechanics simulation. The ternary spin domain {-1, 0, +1} matches the ecosystem's universal ternary representation, enabling direct mapping between Ising energy landscapes and ternary neural network loss surfaces.

## References

- **Original Ising Model**: Ising, E., "Beitrag zur Theorie des Ferromagnetismus," Zeitschrift für Physik, 31, 253-258, 1925.
- **Exact 2D Solution**: Onsager, L., "Crystal Statistics. I. A Two-Dimensional Model with an Order-Disorder Transition," Physical Review, 65(3-4), 117-149, 1944.
- **Blume-Capel Model**: Blume, M., "Theory of the First-Order Magnetic Phase Change in UO₂," Physical Review, 141(2), 517-524, 1966. Capel, H.W., "On the possibility of first-order phase transitions in Ising systems," Physica, 32(5), 966-988, 1966.
- **Metropolis Algorithm**: Metropolis, N. et al., "Equation of State Calculations by Fast Computing Machines," Journal of Chemical Physics, 21(6), 1087-1092, 1953.
- **Monte Carlo Methods**: Newman, M.E.J. & Barkema, G.T., "Monte Carlo Methods in Statistical Physics," Oxford University Press, 1999.
- **Topological Insulators**: Hasan, M.Z. & Kane, C.L., "Colloquium: Topological insulators," Reviews of Modern Physics, 82(4), 3045-3067, 2010.

## License

MIT
