# Projected ED Design — 2026-04-01

## Goal

Implement band-projected exact diagonalization (projected ED): project the two-band model onto the lowest Chern band and perform ED entirely within that reduced Hilbert space. This is a new git branch `projected-ed`; original code is untouched.

## Directory Structure

```
projected_ed/
├── shared_proj/
│   ├── band.jl              # single-particle band energies + Bloch eigenstates
│   ├── proj_interaction.jl  # projected two-body matrix elements W(k1,k2,q)
│   ├── proj_basis.jl        # k-space Fock states + fermion sign + momentum sectors
│   └── proj_hamiltonian.jl  # build H matrix per momentum sector, full diag / Lanczos
└── case_proj/
    └── main.jl              # top-level script
```

Reuse `shared/lattice.jl` and `shared/hoppings.jl` (via `include("../shared/...")`) — no copies.

## Physics

### Step 1 — Single-particle band (band.jl)

For each k-point `k_m` (m = 0..Nk-1), diagonalize the existing 2×2 `get_Hk(k, t1, t3)`:

- `eps[m]`: lower-band eigenvalue ε_m^−
- `uvec[m]`: lower-band Bloch eigenvector [u_{m,A}, u_{m,B}], normalized

Gauge fix: rotate so that `uvec[m][1]` is real positive.

### Step 2 — Projected interaction (proj_interaction.jl)

Form factor for sublattice α:
```
F_α(k', k) = conj(uvec[k'][α]) * uvec[k][α]
```

**V1V2V3 version** — real-space sublattice-resolved Fourier transform:
```
V_αβ(q) = Σ_{ΔR} V(|ΔR + τ_α − τ_β|) × exp(−i·q·(ΔR + τ_α − τ_β))
```
where V(r) = V1, V2, V3 by neighbor distance; τ_A=(0,0), τ_B=(−0.5, √3/2).

**Gate Coulomb version** — dual-gate screened potential:
```
V(q) = (2π / |q|) × tanh(|q| × d),    V(q=0) = 0
```
Parameter `d` is the gate distance (in units of lattice constant). The q=0 term is set to zero (neutralizing background).

Projected two-body matrix element:
```
W(m1, m2, q) = Σ_{α,β} V_αβ(q) × F_α(m1+q, m1) × F_β(m2−q, m2)
```
Precomputed and stored as a 3D array `W[m1, m2, q_label]` (Nk³ entries, Nk=15 → 3375).

### Step 3 — k-space Fock states (proj_basis.jl)

- State: integer bitmask of length Nk (bit m = 1 means k_m occupied)
- Total momentum: `M = (Σ_{occupied m} m) mod Nk`
- Fermion sign for annihilating k_m: `(−1)^(number of occupied bits below m)`
- Enumerate all C(Nk, Np) states, group by M → Nk momentum sectors

### Step 4 — Hamiltonian and diagonalization (proj_hamiltonian.jl)

For each momentum sector M, build the dim×dim matrix:
```
H[i,i] += Σ_{occupied k in state_i} eps[k]          # kinetic diagonal
H[i,j] += (1/Nk) × W(m1,m2,q) × fermion_sign        # interaction off-diagonal
```
Scattering: k1 → k1+q, k2 → k2-q (with k1+q and k2-q in-bounds and satisfying Pauli).

Solve with `eigen(Hermitian(H))` (full diag; upgrade to Lanczos for larger systems).

## Validation

1. **Pure kinetic limit** (V=0 or d→0): ground-state energy = Σ of Np lowest ε_k, error < 1e-10
2. **Trace check**: tr(H) matches sum of all diagonal elements analytically
3. **Particle number check**: ⟨Ψ|Σ_k n_k|Ψ⟩ = Np for all eigenstates
4. **Cross-check with full ED**: on RectLat3x6 (Nk=9, Np=3), projected ED energy ≥ full two-band ED energy (projection is a variational restriction)

## Output Format

Files in `case_proj/`:
```
spectrum_proj_V123.dat    # columns: m  E-E0
spectrum_proj_gate.dat    # columns: m  E-E0
```
Same format as existing `spectrum_Np*.dat` for direct reuse of `plot.jl`.

## Git Strategy

- Branch: `projected-ed` (off `master`)
- Only `projected_ed/` is new; all existing files untouched
- `docs/` directory added to both branches
