# ============================================================
# Projected ED — lowest Chern band
# Lattice: RectLat3x6 (Nk=9) for testing; swap for TiltedLat30 (Nk=15)
# ============================================================
include("../../shared/lattice.jl")
include("../../shared/hoppings.jl")
include("../shared_proj/band.jl")
include("../shared_proj/proj_basis.jl")
include("../shared_proj/proj_interaction.jl")
include("../shared_proj/proj_hamiltonian.jl")

using Printf, LinearAlgebra, Dates

# ── Parameters ──
t1     = 1.0
t3     = 0.2
V1     = 1.0;  V2 = 0.0;  V3 = 0.0
d_gate = 1.0
Np     = 3

println("="^60)
println("Projected ED — lowest Chern band")
println("t1=$t1  t3=$t3  Np=$Np")
println("="^60)

# ── 1. Lattice ──
lat = RectLat3x6()   # swap to TiltedLat30() for production
@printf("Lattice: Ns=%d  Nk=%d\n", lat.Ns, lat.Nuc)
Nk = lat.Nuc

# ── 2. Band structure ──
print("Computing band structure... ")
eps, uvec = compute_band(lat, t1, t3)
@printf("done. Band range: [%.4f, %.4f]\n", minimum(eps), maximum(eps))

# ── 3. Basis ──
print("Building k-space Fock basis C($Nk,$Np)=$(binomial(Nk,Np))... ")
sectors = build_proj_sectors(Nk, Np)
println("done. Sector sizes: ", [length(sectors[M]) for M in 0:Nk-1])

# ── 4. Momentum addition table (shared by both interaction types) ──
print("Building k-addition table... ")
kadd = build_k_add(lat)
println("done.")

# ── 5a. V1V2V3 version ──
println("\n--- V1V2V3 version (V1=$V1, V2=$V2, V3=$V3) ---")
W_v123 = build_Vq_V123(lat, uvec, V1, V2, V3)
ev_v123, _ = solve_proj_ed(sectors, eps, W_v123, kadd, Nk)
E0_v123 = ev_v123[1][2]
open("spectrum_proj_V123.dat", "w") do f
    println(f, "# m  E-E0  [proj ED Nk=$Nk Np=$Np V1=$V1 V2=$V2 V3=$V3 t3=$t3]")
    for (m, e) in ev_v123
        @printf(f, "%d  %.10f\n", m, e - E0_v123)
    end
end
println("Ground state: M=$(ev_v123[1][1])  E0=$(round(E0_v123, digits=6))")
println("Saved: spectrum_proj_V123.dat")

# ── 5b. Gate Coulomb version ──
println("\n--- Gate Coulomb version (d=$d_gate) ---")
W_gate = build_Vq_gate(lat, uvec, d_gate)
ev_gate, _ = solve_proj_ed(sectors, eps, W_gate, kadd, Nk)
E0_gate = ev_gate[1][2]
open("spectrum_proj_gate.dat", "w") do f
    println(f, "# m  E-E0  [proj ED Nk=$Nk Np=$Np d=$d_gate t3=$t3]")
    for (m, e) in ev_gate
        @printf(f, "%d  %.10f\n", m, e - E0_gate)
    end
end
println("Ground state: M=$(ev_gate[1][1])  E0=$(round(E0_gate, digits=6))")
println("Saved: spectrum_proj_gate.dat")

println("\nDone: ", now())
