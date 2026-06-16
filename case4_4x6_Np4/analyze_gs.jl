# ============================================================
# Case 4 基态分析：实空间密度 + BZ 结构因子
#
# 用法：julia analyze_gs.jl
# 输出：
#   rspace_Np4.dat   — 每格点 ⟨n_i⟩
#   sq_fine_Np4.dat  — 精细 BZ 网格 N(q)
#   rspace_Np4.pdf   — 实空间密度图
#   sq_bz_Np4.pdf    — BZ 中 N(q) 热图
# ============================================================
include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/hamiltonian.jl")
include("../shared/solver.jl")
include("../shared/structure_factor.jl")
include("../shared/wavefunction.jl")

using Plots, Printf

# ── 参数 ──
Np = 4; t1 = 1.0; t3 = 0.2; V1 = 1.0; V2 = 0.0; V3 = 0.0

println("="^60)
println("基态分析 — 4×6 sites, Np=$Np")
println("t=$t1  t'=$t3  V1=$V1")
println("="^60)

# ── 构建晶格与基矢 ──
print("构建晶格与基矢... "); flush(stdout)
lat   = RectLat4x6()
basis = gen_basis(lat.Ns, Np)
secs  = build_all_ksectors(basis, lat)
hops  = build_hops(lat, t1, t3, 0.0)
println("done")

# ── 求解全谱（nev=1），确定基态扇区 ──
println("Lanczos (nev=1) 全扇区...")
all_ev = compute_spectrum(secs, lat, hops, V1, V2, V3; nev=1)
E0    = all_ev[1][2]
gs_k  = all_ev[1][1]
n_gs  = lat.Nuc ÷ Np    # 预期简并度 = 3

println("\n最低 $(min(6,length(all_ev))) 个能量（相对基态）:")
for (i, (m, e)) in enumerate(all_ev[1:min(6, end)])
    @printf("  %d. k=%2d  ΔE=%.8f\n", i, m, e - E0)
end
@printf("基态扇区: k=%d,  E0=%.8f\n", gs_k, E0)

# ── 提取基态波函数 ──
print("\n提取基态波函数... "); flush(stdout)
gs_sec = secs[findfirst(s -> s.m == gs_k, secs)]
_, gs_vecs = solve_sector_with_vec(gs_sec, lat, hops, V1, V2, V3; nev=1)
coeffs = expand_to_fock(gs_sec, gs_vecs[1])
println("done  (Fock 空间展开: $(length(coeffs)) 个态)")

# ── 实空间占据数 ──
occ = site_occupation(coeffs, lat)
print_rspace_report(occ, lat, Np)

# 保存
open("rspace_Np$(Np).dat", "w") do f
    println(f, "# site  ix  iy  rx  ry  <n_i>  [4x6 Np=$Np V1=$V1 t'=$t3]")
    for (s, (ix, iy)) in enumerate(lat.sites)
        rx, ry = site_phys_pos(lat, ix, iy)
        @printf(f, "%d  %d  %d  %.6f  %.6f  %.8f\n", s, ix, iy, rx, ry, occ[s])
    end
end
println("实空间密度 → rspace_Np$(Np).dat")

# ── 团簇 k 点结构因子 ──
results_ord, results_srt, norm2, Np_check = structure_factor_survey(gs_sec, gs_vecs[1], lat, Np)
print_sq_report(results_srt, norm2, Np_check, Np, lat.Ns, gs_k)

# ── n(q) = n(-q) 验证 ──
n̄ = Np / lat.Ns
println("\n" * "="^60)
println("验证 n(q) = n(-q)（对所有 12 个团簇 k 点的物理向量）")
println("-"^60)
@printf("  %-4s  %-20s  %-10s  %-10s  %-8s\n",
        "ik", "physical q", "n(+q)", "n(-q)", "差值")
println("-"^60)
for (m, kphys) in zip(lat.ktab, lat.kpoints)
    qx, qy = kphys[1], kphys[2]
    n_pos = structure_factor_at(coeffs, lat,  qx,  qy, n̄)
    n_neg = structure_factor_at(coeffs, lat, -qx, -qy, n̄)
    @printf("  k=%-2d  (%.4f, %.4f)  %10.6f  %10.6f  %+.2e\n",
            m, qx, qy, n_pos, n_neg, n_neg - n_pos)
end
println("="^60)

# ── 精细 BZ 网格 N(q) ──
N1, N2 = 24, 36
@printf("\n精细 BZ N(q) 扫描 (%d×%d = %d 点)...\n", N1, N2, N1*N2)
flush(stdout)
qxs, qys, Sqs = structure_factor_bz(coeffs, lat, N1, N2, n̄)

open("sq_fine_Np$(Np).dat", "w") do f
    println(f, "# qx  qy  N(q)  [BZ grid $(N1)x$(N2), Np=$Np]")
    for i in eachindex(qxs)
        @printf(f, "%.6f  %.6f  %.8f\n", qxs[i], qys[i], Sqs[i])
    end
end
println("精细结构因子 → sq_fine_Np$(Np).dat")

@printf("N(q) 最大值: %.4f at (%.4f, %.4f)\n",
    maximum(Sqs), qxs[argmax(Sqs)], qys[argmax(Sqs)])

# ============================================================
# 绘图
# ============================================================

# ── 图 1：实空间密度图 ──
rxs = [site_phys_pos(lat, ix, iy)[1] for (ix, iy) in lat.sites]
rys = [site_phys_pos(lat, ix, iy)[2] for (ix, iy) in lat.sites]
n̄   = Np / lat.Ns

p1 = scatter(rxs, rys;
    marker_z     = occ,
    color        = cgrad(:RdBu, rev=true),
    markersize   = 16,
    markerstrokewidth = 0.8,
    markerstrokecolor = :black,
    clims        = (0.0, max(maximum(occ)*1.05, 2*n̄)),
    colorbar_title = "⟨n_i⟩",
    xlabel       = "x (a₁)",
    ylabel       = "y",
    title        = "4×6 Np=$Np: 实空间密度 ⟨n_i⟩  (n̄=$(round(n̄,digits=4)))\n" *
                   "t=1, t'=0.2, V₁=1",
    legend       = false,
    aspect_ratio = :equal,
    framestyle   = :box,
    size         = (520, 520),
    left_margin  = 8Plots.mm,
    bottom_margin = 10Plots.mm
)

savefig(p1, "rspace_Np$(Np).pdf")
println("\n实空间密度图 → rspace_Np$(Np).pdf")

# ── 图 2：BZ 中 N(q) 散点图（精细网格 + 团簇 k 点叠加） ──
ck_qxs = [r[3] for r in results_ord]
ck_qys = [r[4] for r in results_ord]
ck_Sqs = [r[2] for r in results_ord]

p2 = scatter(qxs, qys;
    marker_z          = Sqs,
    color             = :hot,
    markersize        = 4,
    markerstrokewidth = 0,
    colorbar_title    = "N(q)",
    xlabel            = "qx",
    ylabel            = "qy",
    title             = "N(q) in BZ — 4×6 Np=$Np\n" *
                        "t=1, t'=0.2, V₁=1  (◇=团簇k点)",
    legend            = false,
    aspect_ratio      = :equal,
    framestyle        = :box,
    size              = (520, 520),
    left_margin       = 8Plots.mm,
    bottom_margin     = 10Plots.mm
)
# 叠加团簇 k 点
scatter!(p2, ck_qxs, ck_qys;
    marker            = :diamond,
    markersize        = 8,
    color             = :white,
    markerstrokecolor = :black,
    markerstrokewidth = 1.2,
    label             = false
)

savefig(p2, "sq_bz_Np$(Np).pdf")
println("BZ 结构因子图 → sq_bz_Np$(Np).pdf")

println("\n完成!")
