# ============================================================
# Case 3：倾斜 (4×4-1) 团簇，30 格点，Np=12（填充 12/30=2/5）
# 参数：t=1, t'=0.2, V1=1, V2=V3=0
#
# 服务器运行：julia --threads 64 main.jl
#   C(30,12) = 86,493,225；15 个扇区每扇区约 5.77M 维
#   主要瓶颈：build_all_ksectors（~186B 次操作）
#             Lanczos 求解（5.77M × krylovdim=200）
# ============================================================
include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/hamiltonian.jl")
include("../shared/solver.jl")
include("../shared/structure_factor.jl")

using Printf, LinearAlgebra, Dates, SparseArrays

# ── BLAS 线程分配：64核 / 15扇区 ≈ 每扇区4线程，提升 dot/axpy 效率 ──
# 扇区数在晶格构建后才知，这里用已知值 15（TiltedLat30 的 Nuc）
let blas_per_sec = max(1, Threads.nthreads() ÷ 15)
    BLAS.set_num_threads(blas_per_sec)
    @printf("BLAS threads/扇区: %d  (Julia threads: %d)\n",
            blas_per_sec, Threads.nthreads())
end

# ── 参数 ──
t1  = 1.0
t3  = 0.2
V1  = 10.0
V2  = 2.0
V3  = 2.0
Np  = 12
nev = 8
kd  = 100    # 降低峰值内存：15×100×5.77M×16B + 31GB CSR ≈ 169 GB（留 ~87 GB 余量）

println("="^60)
println("ED — 倾斜 (4×4-1) 团簇，30 sites，Np=$Np（填充 $(Np)/30 = 2/5）")
println("t=$t1  t'=$t3  V1=$V1  V2=$V2  V3=$V3")
@printf("Julia threads: %d  BLAS threads: %d\n",
        Threads.nthreads(), BLAS.get_num_threads())
println("开始时间: ", now())
println("="^60)

t_total = time()

# ── 1. 构建晶格 ──
print("[1/5] 构建晶格... "); flush(stdout)
t = @elapsed lat = TiltedLat30()
@printf("done  Ns=%d  Nuc=%d  (%.2f s)\n", lat.Ns, lat.Nuc, t)

# ── 2. 生成基矢 ──
@printf("[2/5] 生成基矢 C(%d,%d)=%d... ", lat.Ns, Np, binomial(lat.Ns, Np))
flush(stdout)
t = @elapsed basis = gen_basis(lat.Ns, Np)
@printf("done  (%.2f s)\n", t)
@printf("      内存估算: %.1f MB\n", length(basis) * 8 / 1e6)

# ── 3. 构建动量扇区（逐扇区计时，最大瓶颈）──
println("[3/5] 构建动量扇区（逐扇区计时）...")
println("      注：每扇区需对 $(length(basis)) 个 Fock 态做 $(lat.Nuc) 次平移，是主要瓶颈")
flush(stdout)
secs = KSector[]
t_secs = @elapsed begin
    for (i, m) in enumerate(lat.ktab)
        local t = @elapsed sec = build_ksector(basis, lat, m)
        push!(secs, sec)
        elapsed_total = time() - t_total
        @printf("  k=%2d (%2d/%d): dim=%6d  本扇区 %5.1f s  累计 %5.1f s\n",
                m, i, lat.Nuc, length(sec.reps), t, elapsed_total)
        flush(stdout)
    end
end
total_reps = sum(length(s.reps) for s in secs)
@printf("  总代表元: %d / %d ✓  扇区构造总耗时: %.1f s\n",
        total_reps, length(basis), t_secs)
@printf("  fock2rep 内存估算: %.1f MB\n",
        sum(length(s.fock2rep) for s in secs) * (8+4+16) / 1e6)

# ── 4. CSR 稀疏矩阵能谱（步骤 A 建矩阵 + 步骤 B Lanczos）──
println("\n[4/5] CSR 稀疏矩阵能谱（krylovdim=$(kd)，$(Threads.nthreads()) 线程）...")
println("      内存预估：CSR ~31 GB + Krylov ~$(round(15*kd*5.77e6*16/1e9,digits=0)) GB")
println("      保留各扇区基态波函数，供结构因子直接使用")
flush(stdout)
hops0 = build_hops(lat, t1, t3, 0.0)
t_spec = @elapsed begin
    all_ev, gs_vecs = compute_spectrum_sparse_with_vecs(
        secs, lat, hops0, V1, V2, V3;
        nev=nev, krylovdim=kd, verbose=true)
end
E0 = all_ev[1][2]
@printf("  能谱总耗时: %.1f s (%.2f min)\n", t_spec, t_spec/60)

println("\n  最低 20 个能量（相对基态）:")
for (i,(m,e)) in enumerate(all_ev[1:min(20,end)])
    @printf("    %2d. k=%2d  ΔE=%.8f\n", i, m, e-E0)
end
n15 = min(15, length(all_ev))
gap = length(all_ev) > 15 ? all_ev[16][2]-all_ev[15][2] : NaN
@printf("  前15态散布: %.8f\n", all_ev[n15][2]-E0)
@printf("  15→16 gap:  %.8f\n", gap)

open("spectrum_Np$(Np).dat","w") do f
    println(f,"# k  E-E0  [30sites Np=$Np V1=$V1 V2=$V2 V3=$V3 t'=$t3 CSR]")
    for (m,e) in all_ev
        @printf(f,"%d  %.10f\n", m, e-E0)
    end
end
println("  能谱保存: spectrum_Np$(Np).dat")

# ── 释放 CSR 矩阵（步骤 5 结构因子不需要，释放 ~31 GB）──
# （CSR 变量在 compute_spectrum_sparse_with_vecs 内部，已离开作用域，GC 可回收）
GC.gc()
@printf("  GC 后内存：已释放 CSR 矩阵\n")
flush(stdout)

# ── 5. 结构因子 ──
# 各扇区基态波函数已由步骤 4 保存在 gs_vecs，无需重跑 Lanczos。
println("\n[5/5] 结构因子 N(k) 计算（全扇区，直接使用已保存波函数）...")
flush(stdout)

gs_k   = all_ev[1][1]
gs_sec = secs[findfirst(s -> s.m == gs_k, secs)]

# ── 5a. 全局基态扇区的 N(q) ──
print("  全局基态 (k=$gs_k) 结构因子... "); flush(stdout)
t_sq = @elapsed begin
    gs_vec = gs_vecs[gs_k]
    results_ord, results_srt, norm2, Np_check =
        structure_factor_survey(gs_sec, gs_vec, lat, Np)
end
@printf("done (%.1f s)\n", t_sq)
print_sq_report(results_srt, norm2, Np_check, Np, lat.Ns, gs_k)

open("sq_Np$(Np).dat","w") do f
    println(f,"# ik  N(k)  kx  ky  [30sites Np=$Np V1=$V1 t'=$t3]")
    println(f,"# 基态扇区 k=$gs_k  E0=$(round(E0,digits=8))")
    for (m, sq, kx, ky) in results_ord
        @printf(f,"%d  %.6f  %.6f  %.6f\n", m, sq, kx, ky)
    end
end
println("  结构因子保存: sq_Np$(Np).dat")

# ── 5b. 所有扇区基态的 N(q)（无额外 Lanczos 开销）──
print("  全扇区结构因子矩阵... "); flush(stdout)
t_all = @elapsed open("sq_all_Np$(Np).dat","w") do f
    println(f,"# sector_k  q  N(q)  kx  ky  [30sites Np=$Np V1=$V1 t'=$t3]")
    println(f,"# 行：扇区 k（本征态来自该扇区基态）  列：动量 q")
    for sec in secs
        haskey(gs_vecs, sec.m) || continue
        res_ord, _, _, _ = structure_factor_survey(sec, gs_vecs[sec.m], lat, Np)
        for (m, sq, kx, ky) in res_ord
            @printf(f,"%d  %d  %.6f  %.6f  %.6f\n", sec.m, m, sq, kx, ky)
        end
    end
end
@printf("done (%.1f s)\n", t_all)
println("  全扇区结构因子保存: sq_all_Np$(Np).dat")

@printf("\n总耗时: %.1f s (%.2f min)\n", time()-t_total, (time()-t_total)/60)
println("完成时间: ", now())
println("完成!")
