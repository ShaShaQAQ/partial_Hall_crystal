# ============================================================
# Case 5：4×6 sites（4×3 原胞），24 格点，Np=4，V1=0（无相互作用）
# 填充：n̄ = 4/24 = 1/6；ν = 4/12 = 1/3（下 Chern 带）
# 参数：t=1, t'=0.2, V1=0, V2=V3=0
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
t1      = 1.0
t3      = 0.2
V1      = 0.0
V2      = 0.0
V3      = 0.0
Np      = 4
nev     = 4

println("="^60)
println("ED — 4×6 sites（4×3 原胞），Ns=24，Np=$Np")
println("n̄=$(Np)/24=$(round(Np/24, digits=4))  ν=$(Np)/12=$(round(Np/12, digits=4))")
println("t=$t1  t'=$t3  V1=$V1  V2=$V2  V3=$V3")
println("="^60)

# ── 构建晶格 ──
print("构建晶格... "); flush(stdout)
lat = RectLat4x6()
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
println("\n零通量能谱（Lanczos，保留各扇区基态波函数）...")
hops0 = build_hops(lat, t1, t3, 0.0)
all_ev, gs_vecs = compute_spectrum_with_vecs(secs, lat, hops0, V1, V2, V3; nev=nev)
E0 = all_ev[1][2]

println("\n最低 $(min(15,length(all_ev))) 个能量（相对基态）:")
for (i,(m,e)) in enumerate(all_ev[1:min(15,end)])
    @printf("  %2d. k=%d  ΔE=%.8f\n", i, m, e-E0)
end

n3 = min(3, length(all_ev))
gap = length(all_ev) > 3 ? all_ev[4][2]-all_ev[3][2] : NaN
@printf("\n前3态散布: %.8f\n", all_ev[n3][2]-E0)
@printf("3→4 gap:   %.8f\n", gap)

open("spectrum_Np$(Np).dat","w") do f
    println(f,"# k  E-E0  [4x6 Np=$Np V1=$V1 V2=$V2 V3=$V3 t'=$t3]")
    for (m,e) in all_ev
        @printf(f,"%d  %.10f\n", m, e-E0)
    end
end
println("能谱保存: spectrum_Np$(Np).dat")

# ── 结构因子 N(k) ──
# 直接使用 Lanczos 阶段已保存的波函数，无需重跑。
println("\n结构因子 N(k) 计算（CDW 诊断）...")
gs_k   = all_ev[1][1]
gs_sec = secs[findfirst(s -> s.m == gs_k, secs)]

# ── 全局基态扇区 ──
gs_vec = gs_vecs[gs_k]
results_ord, results_srt, norm2, Np_check =
    structure_factor_survey(gs_sec, gs_vec, lat, Np)
print_sq_report(results_srt, norm2, Np_check, Np, lat.Ns, gs_k)

open("sq_Np$(Np).dat","w") do f
    println(f,"# ik  N(k)  kx  ky  [4x6 Np=$Np V1=$V1 t'=$t3]")
    println(f,"# 基态扇区 k=$gs_k  E0=$(round(E0,digits=8))")
    for (m, sq, kx, ky) in results_ord
        @printf(f,"%d  %.6f  %.6f  %.6f\n", m, sq, kx, ky)
    end
end
println("结构因子保存: sq_Np$(Np).dat")

# ── 全扇区结构因子矩阵 ──
open("sq_all_Np$(Np).dat","w") do f
    println(f,"# sector_k  q  N(q)  kx  ky  [4x6 Np=$Np V1=$V1 t'=$t3]")
    for sec in secs
        haskey(gs_vecs, sec.m) || continue
        res_ord, _, _, _ = structure_factor_survey(sec, gs_vecs[sec.m], lat, Np)
        for (m, sq, kx, ky) in res_ord
            @printf(f,"%d  %d  %.6f  %.6f  %.6f\n", sec.m, m, sq, kx, ky)
        end
    end
end
println("全扇区结构因子保存: sq_all_Np$(Np).dat")

println("\n完成!")
