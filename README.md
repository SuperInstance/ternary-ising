# ternary-ising

**The ternary Ising model — where 0 is a topological insulator.**

The Ising model is *the* workhorse of statistical mechanics: spins on a lattice, interacting with their neighbors, ordered at low temperature and chaotic at high temperature. The classic binary Ising model (up/down) has a famous phase transition at a critical temperature.

This crate implements a *ternary* Ising model where spins take values `{-1, 0, +1}`. The 0 spin is special — it's a **topological insulator**. It doesn't participate in magnetic ordering, doesn't feel the coupling force, and screens interactions between +1 and -1 spins. The result: **no phase transition**. The ternary Ising model stays disordered at every temperature. The 0 state kills long-range order.

This was discovered in the spiral experiments and documented as a fundamental property of ternary systems.

## What's Inside

- **`TernaryIsing`** — 2D lattice with spins, temperature (β), and coupling (J)
- **`Spin::Down / Zero / Up`** — the three spin states
- **`local_energy(r, c)`** — energy at a site from neighbor interactions
- **`total_energy()`** — system energy
- **`magnetization()`** — net alignment (always near 0 for ternary!)
- **`metropolis_step(rng)`** — single Monte Carlo update step
- **`simulate(steps, rng)`** — full simulation with energy/magnetization history
- **`susceptibility(history)`** — magnetic susceptibility from magnetization fluctuations

## Quick Example

```rust
use ternary_ising::*;

// 20x20 lattice, moderate temperature
let mut model = TernaryIsing::new(20, 20, 1.0, 1.0);

// Run Monte Carlo simulation
let mut rng = || 0.5; // your RNG
model.simulate(1000, &mut rng);

// Key finding: magnetization ≈ 0 regardless of temperature
let m = model.magnetization();
println!("Magnetization: {:.3}", m); // ~0.05 — no ordering!

// Compare: binary Ising would show M → ±1 at low temperature
// Ternary: the 0 insulator prevents long-range order
```

## The Insight

**The 0 state screens everything.** In binary Ising, spins align or anti-align. In ternary, the 0 spin is neutral — it breaks the chain of interactions. Energy can't propagate through 0 cells. This means no phase transition, no spontaneous symmetry breaking, no ordered state. The ternary Ising model lives in a *single phase* — perpetually disordered.

This is the same mechanism that makes ternary flocking fail, ternary Kuramoto fail, and ternary systems resist monoculture. The 0 is a universal screen.

**Use cases:**
- **Statistical mechanics research** — study the effect of neutral states on phase transitions
- **Material science** — model ternary alloys and topological insulators
- **Agent systems** — understand why ternary populations resist ordering
- **Education** — the simplest non-trivial Ising model
- **Complex systems** — how neutral particles affect collective behavior

## See Also

- **ternary-life** — ternary Game of Life (another lattice dynamics model)
- **ternary-kuramoto** — ternary synchronization (also fails due to 0)
- **ternary-fire** — forest fire model (ternary states, different dynamics)
- **ternary-sandpile** — self-organized criticality on ternary grids
- **ternary-percolation** — percolation through ternary lattices

## Install

```bash
cargo add ternary-ising
```

## License

MIT
