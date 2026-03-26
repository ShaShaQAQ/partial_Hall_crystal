# ============================================================
# Case 1：3×6 矩形超胞，18 格点，Np=6（填充 1/3）
# 参数：t=1, t'=0.2, V1=1, V2=V3=0
# 对应论文 Fig.1(a) 相图（FCI vs CDW）
# ============================================================
push!(LOAD_PATH, @__DIR__)
include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/hamiltonian.jl")
include("../shared/solver.jl")
include("../shared/structure_factor.jl")

using Printf

# ── 参数 ──
t1      = 1.0
t3      = 0.2
V1      = 1.0
V2      = 0.0
V3      = 0.0
Np      = 6
nev     = 4
n_phi   = 60
phi_max = 6π      # 0→3 个通量量子（FCI ν=2/3 周期 = 3×2π）

println("="^60)
println("ED — 3×6 矩形超胞，18 sites，Np=$Np（填充 $(Np)/18 = 1/3）")
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

# ── 零通量能谱 ──
println("\n零通量能谱（Lanczos）...")
hops0 = build_hops(lat, t1, t3, 0.0)
all_ev = compute_spectrum(secs, lat, hops0, V1, V2, V3; nev=nev)
E0 = all_ev[1][2]

println("\n最低 $(min(12,length(all_ev))) 个能量（相对基态）:")
for (i,(m,e)) in enumerate(all_ev[1:min(12,end)])
    @printf("  %2d. k=%d  ΔE=%.8f\n", i, m, e-E0)
end

gap = length(all_ev) > 4 ? all_ev[4][2]-all_ev[3][2] : NaN
@printf("\n3→4 gap: %.8f\n", gap)

# 保存能谱
open("spectrum_Np$(Np).dat","w") do f
    println(f,"# k  E-E0  [3x6 Np=$Np V1=$V1 V2=$V2 V3=$V3 t'=$t3]")
    for (m,e) in all_ev
        @printf(f,"%d  %.10f\n", m, e-E0)
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

    open("sq_Np$(Np).dat","w") do f
        println(f,"# ik  N(k)  kx  ky  [3x6 Np=$Np V1=$V1 t'=$t3]")
        println(f,"# 基态扇区 k=$gs_k  E0=$(round(E0,digits=8))")
        println(f,"# CDW 峰预期在 ik=5,7 (±K)")
        for (m, sq, kx, ky) in results_ord
            @printf(f,"%d  %.6f  %.6f  %.6f\n", m, sq, kx, ky)
        end
    end
    println("结构因子保存: sq_Np$(Np).dat")
end

# ── 谱流 ──
println("\n谱流扫描 φ: 0 → $(round(phi_max/π,digits=1))π ($n_phi 点)...")
phis, all_ens, E0f = spectral_flow(secs, lat, t1, t3, V1, V2, V3;
                                    nev=nev, n_phi=n_phi, phi_max=phi_max)

open("flux_Np$(Np).dat","w") do f
    println(f,"# phi/2pi  E1-E0  E2-E0 ... (nev*Nuc 列)")
    for (i,phi) in enumerate(phis)
        @printf(f,"%.6f", phi/(2π))
        for e in all_ens[i][1:min(lat.Nuc, length(all_ens[i]))]
            @printf(f,"  %.8f", e-E0f)
        end
        println(f)
    end
end
println("谱流保存: flux_Np$(Np).dat")
println("\n完成!")
