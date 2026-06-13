# Ternary Ising — Ising Model with Three-State Spins and Topological Insulators

**Ternary Ising** simulates the Ising model on a 2D lattice where each site carries a spin in {-1, 0, +1}. The zero spin acts as a **topological insulator**: it doesn't interact with neighbors, breaking interaction chains and creating novel phase behavior distinct from the binary Ising model. The crate provides Metropolis-Hastings dynamics, energy computation, and critical temperature detection.

## Why It Matters

The binary Ising model is the most studied model in statistical mechanics, but real systems often have three states: spin-up, spin-down, and spin-zero (disordered). Adding the zero state fundamentally changes the physics: it introduces a **dilution transition** where interaction chains break, creating percolation-like phenomena alongside the magnetic transition. For ternary neural networks, this is directly relevant — zero-weight connections don't propagate signals, and understanding the percolation threshold determines when a network becomes too sparse to function.

## How It Works

### Hamiltonian

The system energy is:

```
E = -J · Σ sᵢ · sⱼ   (nearest-neighbor sum)
```

where sᵢ ∈ {-1, 0, +1}. When sᵢ = 0, the bond contributes zero energy regardless of neighbors — the site is "dark" to interactions.

### Local Energy

`local_energy(r, c) = -J · s(r,c) · Σ s(neighbors)`. O(1) per site with 4-connected neighbors and periodic boundary conditions.

### Metropolis-Hastings

Each Monte Carlo step:
1. Select random site (r, c)
2. Propose new spin s' ∈ {-1, 0, +1} \ {current}
3. Compute ΔE = E(new) - E(old)
4. Accept with probability min(1, exp(-β·ΔE)) where β = 1/(kT)

At high temperature (β → 0): all moves accepted, system is disordered.
At low temperature (β → ∞): only energy-lowering moves accepted, system orders.
At critical βc: scale-free fluctuations, universality.

### Magnetization and Susceptibility

```
M = |⟨s⟩|                    (order parameter)
χ = (⟨M²⟩ - ⟨M⟩²) / T       (susceptibility, peaks at Tc)
```

### Periodic Boundaries

The lattice wraps toroidally: site (0, y) neighbors (width-1, y), and (x, 0) neighbors (x, height-1). This eliminates boundary effects and ensures translational invariance.

### Zero-State Physics

The key novelty: when a site has spin 0, its bonds are "cut." If enough sites are 0, the lattice percolates — clusters of interacting spins become disconnected. The percolation threshold for site dilution on a square lattice is ≈ 0.407.

## Quick Start

```rust
use ternary_ising::{TernaryIsing, Spin};

let mut model = TernaryIsing::new(20, 20, 0.5, 1.0); // beta=0.5, J=1.0

// Run Monte Carlo steps
for _ in 0..10000 {
    let r = (rand::random::<usize>()) % 20;
    let c = (rand::random::<usize>()) % 20;
    // ... propose, compute ΔE, accept/reject
}

println!("Energy: {:.2}", model.total_energy());
```

```bash
cargo add ternary-ising
```

## API

| Type / Function | Description |
|---|---|
| `Spin` | `Down(-1)`, `Zero(0)`, `Up(+1)` with `is_insulator()` |
| `TernaryIsing` | 2D lattice: `new(w, h, beta, j)`, `local_energy()`, `total_energy()` |
| `neighbors(r, c)` | 4-connected with periodic boundaries |

## Architecture Notes

The Ising model is the physics foundation of **SuperInstance** agent dynamics. Agents in state 0 (analogous to spin 0) don't participate in fleet interactions — they're insulators that break coordination chains. The γ + η = C conservation maps to the Ising energy: γ = ordered spins (low energy), η = disordered spins (high energy). See [Architecture](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## References

- Ising, Ernst. "Beitrag zur Theorie des Ferromagnetismus," *Z. Physik*, 31, 1925.
- Baxter, Rodney. *Exactly Solved Models in Statistical Mechanics*, Academic Press, 1982.
- Stauffer, Dietrich & Aharony, Amnon. *Introduction to Percolation Theory*, 2nd ed., CRC, 1994.

## License

MIT
