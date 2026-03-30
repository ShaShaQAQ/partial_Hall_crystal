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
kd  = 60

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
@printf("done  Ns=%d  Nuc=%d  (%.2f s)  RSS=%.2f GB\n", lat.Ns, lat.Nuc, t, mem_rss_gb())

# ── 2. 生成基矢 ──
@printf("[2/5] 生成基矢 C(%d,%d)=%d... ", lat.Ns, Np, binomial(lat.Ns, Np))
flush(stdout)
t = @elapsed basis = gen_basis(lat.Ns, Np)
@printf("done  (%.2f s)\n", t)
@printf("[MEM] 生成基矢后 RSS = %.2f GB\n", mem_rss_gb())
flush(stdout)

# ── 3. 流水线：逐扇区 ksector → CSR → 立即释放 orbit_data/fock2rep ──
# 内存优化关键：不再批量构建所有 ksector，而是每扇区建完 CSR 立即释放大内存。
# 峰值从 ~39 GB（旧）降至 ~11 GB/扇区（新）。
println("[3/5] 流水线构建扇区+CSR（扇区 $(seg_start)–$(seg_end)）...")
println("      每扇区：build_ksector → build_sparse_H → empty!(orbit_data/fock2rep) → GC")
flush(stdout)
t = @elapsed hops0 = build_hops(lat, t1, t3, 0.0)
@printf("  build_hops %.2f s  RSS=%.2f GB\n", t, mem_rss_gb())
H_csrs    = SparseMatrixCSC{ComplexF64,Int32}[]
secs_local = KSector[]
t_pipeline = @elapsed begin
    for (i, m) in enumerate(lat.ktab)
        i-1 < seg_start && continue
        i-1 > seg_end   && continue

        # 建 ksector
        t_ks = @elapsed sec = build_ksector(basis, lat, m)
        @printf("  k=%2d  dim=%d  ksector %.1f s  RSS=%.2f GB\n",
                m, length(sec.reps), t_ks, mem_rss_gb())
        flush(stdout)

        # 建 CSR（使用 orbit_data + fock2rep）
        t_csr = @elapsed H = build_sparse_H(sec, lat, hops0, V1, V2, V3)
        @printf("  k=%2d  CSR %.1f s  nnz=%d  RSS=%.2f GB\n",
                m, t_csr, nnz(H), mem_rss_gb())
        flush(stdout)

        # 立即释放不再需要的大内存
        empty!(sec.fock2rep)
        empty!(sec.orbit_data)
        GC.gc()
        @printf("  k=%2d  GC 后 RSS=%.2f GB\n", m, mem_rss_gb())
        flush(stdout)

        push!(H_csrs, H)
        push!(secs_local, sec)
    end
end

# 释放 basis（Lanczos 不再需要）
empty!(basis); GC.gc()
total_reps = sum(length(s.reps) for s in secs_local)
@printf("  流水线完成  %d 扇区  代表元合计: %d  耗时: %.1f s\n",
        length(secs_local), total_reps, t_pipeline)
@printf("[MEM] basis 释放后 RSS = %.2f GB\n", mem_rss_gb())
flush(stdout)

# ── 4. 并行 Lanczos（对预建 CSR 矩阵）──
println("\n[4/5] 并行 Lanczos（krylovdim=$(kd)，$(Threads.nthreads()) 线程）...")
flush(stdout)
t_spec = @elapsed begin
    all_ev, gs_vecs = lanczos_sparse_sectors(
        secs_local, H_csrs;
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
@printf("[MEM] 能谱写文件后 RSS = %.2f GB\n", mem_rss_gb())
flush(stdout)

GC.gc()
@printf("[MEM] GC 后 RSS = %.2f GB\n", mem_rss_gb())
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
