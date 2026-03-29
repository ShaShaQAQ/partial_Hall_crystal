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
using JLD2

# ── 多节点参数：本节点负责的扇区索引范围（0-based，闭区间）──
seg_start = parse(Int, get(ENV, "SECTOR_START", "0"))
seg_end   = parse(Int, get(ENV, "SECTOR_END",   "14"))
seg_id    = seg_start   # 用于命名输出文件

# ── BLAS 线程分配：按本节点实际扇区数分配，提升 dot/axpy 效率 ──
let n_local = seg_end - seg_start + 1,
    blas_per_sec = max(1, Threads.nthreads() ÷ n_local)
    BLAS.set_num_threads(blas_per_sec)
    @printf("本节点扇区数: %d  BLAS threads/扇区: %d  (Julia threads: %d)\n",
            n_local, blas_per_sec, Threads.nthreads())
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
@printf("[MEM] 启动后基线 RSS = %.2f GB\n", mem_rss_gb())
flush(stdout)

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
@printf("[MEM] 生成基矢后 RSS = %.2f GB\n", mem_rss_gb())
flush(stdout)

# ── 3. 构建动量扇区（只建本节点负责的扇区）──
println("[3/5] 构建动量扇区（只建扇区 $(seg_start)–$(seg_end)，跳过其余）...")
println("      注：每扇区需对 $(length(basis)) 个 Fock 态做 $(lat.Nuc) 次平移，是主要瓶颈")
flush(stdout)
secs_local = KSector[]
t_secs = @elapsed begin
    for (i, m) in enumerate(lat.ktab)
        i-1 < seg_start && continue
        i-1 > seg_end   && continue
        local t = @elapsed sec = build_ksector(basis, lat, m)
        push!(secs_local, sec)
        elapsed_total = time() - t_total
        @printf("  k=%2d (%2d/%d): dim=%6d  本扇区 %5.1f s  累计 %5.1f s\n",
                m, i, lat.Nuc, length(sec.reps), t, elapsed_total)
        flush(stdout)
    end
end
total_reps = sum(length(s.reps) for s in secs_local)
@printf("  本节点扇区 %d–%d  代表元: %d  构造耗时: %.1f s\n",
        seg_start, seg_end, total_reps, t_secs)
@printf("  fock2rep 内存估算: %.1f MB（将在 CSR 建完后释放）\n",
        sum(length(s.fock2rep) for s in secs_local) * (8+4+16) / 1e6)
@printf("  orbit_data 内存估算: %.1f GB（轨道数×平均轨道长度×24B）\n",
        sum(sum(length(od) for od in s.orbit_data) for s in secs_local) * 24 / 1e9)
@printf("[MEM] 扇区构建完成后 RSS = %.2f GB\n", mem_rss_gb())
flush(stdout)

# ── 4. CSR 稀疏矩阵能谱（步骤 A 建矩阵 + 步骤 B Lanczos）──
println("\n[4/5] CSR 稀疏矩阵能谱（krylovdim=$(kd)，$(Threads.nthreads()) 线程）...")
let n_loc = length(secs_local)
    @printf("      内存预估：CSR ~%.0f GB + Krylov ~%.0f GB\n",
            n_loc*2.1, n_loc*kd*5.77e6*16/1e9)
end
println("      保留各扇区基态波函数，供结构因子直接使用")
flush(stdout)
hops0 = build_hops(lat, t1, t3, 0.0)
t_spec = @elapsed begin
    all_ev, gs_vecs = compute_spectrum_sparse_with_vecs(
        secs_local, lat, hops0, V1, V2, V3;
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
@printf("  GC 后 RSS = %.2f GB（已释放 CSR 矩阵+Krylov向量）\n", mem_rss_gb())
flush(stdout)

# ── 5. 保存中间结果（供 merge.jl 使用）──
mkpath("output")
out_file = "output/partial_$(seg_id).jld2"
jldsave(out_file; ev_pairs=all_ev, gs_vecs=gs_vecs,
        Np=Np, V1=V1, V2=V2, V3=V3, t1=t1, t3=t3)
@printf("  中间结果保存: %s\n", out_file)
@printf("[MEM] JLD2 保存后 RSS = %.2f GB\n", mem_rss_gb())

@printf("\n总耗时: %.1f s (%.2f min)\n", time()-t_total, (time()-t_total)/60)
println("完成时间: ", now())
println("完成!")
