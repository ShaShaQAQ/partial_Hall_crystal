# ============================================================
# Case 1：3×6 矩形超胞，18 格点，Np=6（填充 6/18 = 1/3）
# 参数：t=1, t'=0.22, V1=1.6, V2=V3=0
# ============================================================
include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/hamiltonian.jl")
include("../shared/solver.jl")
include("../shared/structure_factor.jl")

using Printf

# ── 参数 ──
t1  = 1.0
t3  = 0.22
V1  = 1.6
V2  = 0.0
V3  = 0.0
Np  = 6
nev = 4

println("="^60)
println("ED — 3×6 超胞，Ns=18，Np=$Np（ν=$(Np÷gcd(Np,9))/$(9÷gcd(Np,9)) 下 Chern 带）")
println("t=$t1  t'=$t3  V1=$V1  V2=$V2  V3=$V3")
println("="^60)

# ── 构建晶格 ──
print("构建晶格... "); flush(stdout)
lat = RectLat3x6()
println("done  Ns=$(lat.Ns)  Nuc=$(lat.Nuc)  k-points=$(length(lat.ktab))")

# ── 生成基矢 ──
@printf("生成基矢 C(%d,%d)=%d... ", lat.Ns, Np, binomial(lat.Ns, Np))
flush(stdout)
basis = gen_basis(lat.Ns, Np)
println("done")

# ── 构建动量扇区 ──
println("构建动量扇区...")
secs = build_all_ksectors(basis, lat)
total_reps = sum(length(s.reps) for s in secs)
println("总代表元: $total_reps / $(length(basis)) ✓")
for s in secs
    @printf("  k=%2d: dim=%d\n", s.m, length(s.reps))
end

# ── 零通量能谱 ──
println("\n零通量能谱（Lanczos）...")
hops0 = build_hops(lat, t1, t3, 0.0)
all_ev = compute_spectrum(secs, lat, hops0, V1, V2, V3; nev=nev)
E0 = all_ev[1][2]

n_gs = lat.Nuc ÷ gcd(Np, lat.Nuc)   # 预期简并度 = 3

println("\n最低 $(min(12, length(all_ev))) 个能量（相对基态）:")
for (i, (m, e)) in enumerate(all_ev[1:min(12, end)])
    @printf("  %2d. k=%d  ΔE=%.8f\n", i, m, e - E0)
end

n3  = min(n_gs, length(all_ev))
gap = length(all_ev) > n_gs ? all_ev[n_gs+1][2] - all_ev[n_gs][2] : NaN
@printf("\n前%d态散布: %.8f\n", n_gs, all_ev[n3][2] - E0)
@printf("%d→%d gap:  %.8f\n", n_gs, n_gs+1, gap)

open("spectrum_Np$(Np).dat", "w") do f
    println(f, "# k  E-E0  [3x6 Np=$Np V1=$V1 V2=$V2 V3=$V3 t'=$t3]")
    for (m, e) in all_ev
        @printf(f, "%d  %.10f\n", m, e - E0)
    end
end
println("能谱保存: spectrum_Np$(Np).dat")

# ── 结构因子 N(k) ──
println("\n结构因子 N(k) 计算（CDW 诊断）...")
gs_k   = all_ev[1][1]
gs_sec = secs[findfirst(s -> s.m == gs_k, secs)]
_, gs_vecs = solve_sector_with_vec(gs_sec, lat, hops0, V1, V2, V3; nev=1)

if !isempty(gs_vecs)
    results_ord, results_srt, norm2, Np_check =
        structure_factor_survey(gs_sec, gs_vecs[1], lat, Np)
    print_sq_report(results_srt, norm2, Np_check, Np, lat.Ns, gs_k)

    open("sq_Np$(Np).dat", "w") do f
        println(f, "# ik  N(k)  kx  ky  [3x6 Np=$Np V1=$V1 t'=$t3]")
        println(f, "# 基态扇区 k=$gs_k  E0=$(round(E0, digits=8))")
        for (m, sq, kx, ky) in results_ord
            @printf(f, "%d  %.6f  %.6f  %.6f\n", m, sq, kx, ky)
        end
    end
    println("结构因子保存: sq_Np$(Np).dat")
end

println("\n完成!")
