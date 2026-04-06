# ============================================================
# validate_4x6.jl — 与 case4/case5 全 ED 交叉验证
# RectLat4x6, Np=4 (ν=1/3), case4: V1=1, case5: V1=0
#
# 物理逻辑：FCI/费米液体是单带物理，因此投影 ED 的低能谱
# 应当与全 2-band ED 吻合（投影能量 ≥ 全 ED 能量）。
# ============================================================
include("../../shared/lattice.jl")
include("../../shared/hoppings.jl")
include("../shared_proj/band.jl")
include("../shared_proj/proj_basis.jl")
include("../shared_proj/proj_interaction.jl")
include("../shared_proj/proj_hamiltonian.jl")

using Printf, LinearAlgebra

t1 = 1.0;  t3 = 0.2;  Np = 4

lat  = RectLat4x6()
Nk   = lat.Nuc   # 12

eps, uvec = compute_band(lat, t1, t3)
sectors   = build_proj_sectors(Nk, Np)
kadd      = build_k_add(lat)

@printf("Lattice: Nk=%d  Np=%d  C(%d,%d)=%d\n", Nk, Np, Nk, Np, binomial(Nk,Np))
@printf("Sector sizes: %s\n", string([length(sectors[m]) for m in 0:Nk-1]))
@printf("Band range: [%.6f, %.6f]\n\n", minimum(eps), maximum(eps))

# ── 读取全 ED 数据（相对基态能量）──
function read_full_ed(path)
    evs = Tuple{Int,Float64}[]
    open(path) do f
        for line in eachline(f)
            startswith(line, "#") && continue
            parts = split(strip(line))
            length(parts) == 2 || continue
            push!(evs, (parse(Int, parts[1]), parse(Float64, parts[2])))
        end
    end
    sort!(evs, by=x->x[2])
    return evs
end

# ── 运行一组参数的投影 ED + 对比 ──
function run_and_compare(V1::Float64, full_ed_path::String, label::String)
    println("="^60)
    println("$label  (V1=$V1, V2=V3=0, t3=$t3, Np=$Np)")
    println("="^60)

    # 投影 ED
    W   = build_Vq_V123(lat, uvec, V1, 0.0, 0.0)
    evp, _ = solve_proj_ed(sectors, eps, W, kadd, Nk)
    E0p = evp[1][2]

    # 全 ED
    evf = read_full_ed(full_ed_path)
    E0f_rel = evf[1][2]   # should be 0 (already relative)
    # the .dat stores E-E0_full; we reconstruct absolute E0_full from kinetic sum check
    # Actually we only need relative energies for comparison

    N_show = min(15, length(evp), length(evf))
    @printf("\n%-5s  %-6s  %-14s  %-14s  %-12s\n",
            "rank", "k_proj", "ΔE_proj", "ΔE_full", "deviation")
    println("-"^60)
    for i in 1:N_show
        m_p, e_p = evp[i]
        m_f, e_f = evf[i]
        dev = e_p - E0p - e_f   # both relative to their own E0
        @printf("%-5d  %-6d  %-14.8f  %-14.8f  %+.2e\n",
                i, m_p, e_p - E0p, e_f, dev)
    end

    # FCI 诊断：1/3 填充应有3重简并基态
    n_degen = count(e -> abs(e[2] - E0p) < 1e-4, evp)   # 容差 1e-4
    spread3 = length(evp) >= 3 ? evp[3][2] - E0p : NaN
    gap34   = length(evp) >= 4 ? evp[4][2] - evp[3][2] : NaN
    @printf("\nProj ED: E0=%.6f  sector M=%d\n", E0p, evp[1][1])
    @printf("  前3态散布: %.8f\n", spread3)
    @printf("  3→4 gap:   %.8f\n", gap34)

    # V1=0：基态应等于4个最低单粒子能量之和
    if V1 == 0.0
        E0_kinetic = sum(sort(eps)[1:Np])
        @printf("  纯动能验证: E0_proj=%.8f  sum4_lowest=%.8f  err=%.2e\n",
                E0p, E0_kinetic, abs(E0p - E0_kinetic))
    end

    # 保存投影 ED 能谱
    outname = "spectrum_proj_4x6_V$(round(Int,V1)).dat"
    open(outname, "w") do f
        println(f, "# k  E-E0  [proj ED RectLat4x6 Np=$Np V1=$V1 t3=$t3]")
        for (m, e) in evp
            @printf(f, "%d  %.10f\n", m, e - E0p)
        end
    end
    println("  Saved: $outname")

    return evp, evf
end

evp4, evf4 = run_and_compare(1.0,
    "../../case4_4x6_Np4/spectrum_Np4.dat",
    "Case 4 (V1=1)")

evp5, evf5 = run_and_compare(0.0,
    "../../case5_4x6_Np4_V0/spectrum_Np4.dat",
    "Case 5 (V1=0)")

println("\n" * "="^60)
println("交叉验证完成。")
println("若投影 ED ≈ 全 ED 低能谱，则 FCI/费米液体为单带物理。✓")
