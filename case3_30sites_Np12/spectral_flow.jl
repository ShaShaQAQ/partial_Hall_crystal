# ============================================================
# 谱流计算：φ_y 从 0 到 4π 扫描（2 个磁通量子），追踪能级演化
#
# 策略：matrix-free Hv!（无 CSR 矩阵），按 k-sector 分节点
#   - 本节点只建 SECTOR_START..SECTOR_END 这几个 KSector
#   - 对全部 N_phi 个 φ 点都跑（所有磁通量）
#   - 每个 φ 点只重建 hops（~1 s），对角元 diag_H 预计算一次
#
# 服务器运行（单节点负责扇区 0-4）：
#   SECTOR_START=0 SECTOR_END=4 julia --threads 64 spectral_flow.jl
# ============================================================
include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/hamiltonian.jl")
include("../shared/solver.jl")

using Printf, LinearAlgebra, Dates, JLD2
using KrylovKit

# ── 节点参数：本节点负责的 k-sector 索引范围（0-based，闭区间）──
seg_start = parse(Int, get(ENV, "SECTOR_START", "0"))
seg_end   = parse(Int, get(ENV, "SECTOR_END",   "4"))

# ── 物理参数 ──
t1  = 1.0
t3  = 0.2
V1  = 10.0
V2  = 2.0
V3  = 2.0
Np  = 12
nev = 2
kd  = 60

# ── φ 网格：[0, 4π) 均匀分 N_phi 个点（扫 2 个磁通量子）──
N_phi    = 40
phi_grid = [2π * 2 * i / N_phi for i in 0:(N_phi-1)]

println("="^60)
println("谱流计算 — 倾斜 (4×4-1) 团簇，30 sites，Np=$Np (ν=2/5)")
@printf("本节点扇区: %d – %d  (共 %d 个)\n",
        seg_start, seg_end, seg_end - seg_start + 1)
@printf("φ 点数: %d  范围: 0 → %.4f (%.1f×2π)\n",
        N_phi, phi_grid[end], phi_grid[end]/(2π))
@printf("t=$t1  t'=$t3  V1=$V1  V2=$V2  V3=$V3\n")
@printf("nev=%d  krylovdim=%d\n", nev, kd)
@printf("Julia threads: %d  BLAS threads: %d\n",
        Threads.nthreads(), BLAS.get_num_threads())
println("开始时间: ", now())
println("="^60)

t_total = time()
@printf("[MEM] 启动后基线 RSS = %.2f GB\n", mem_rss_gb()); flush(stdout)

# ── 1. 构建晶格 ──
print("[1/4] 构建晶格... "); flush(stdout)
t = @elapsed lat = TiltedLat30()
@printf("done  Ns=%d  Nuc=%d  (%.2f s)  RSS=%.2f GB\n",
        lat.Ns, lat.Nuc, t, mem_rss_gb())

# ── 2. 生成基矢 ──
@printf("[2/4] 生成基矢 C(%d,%d)=%d... ", lat.Ns, Np, binomial(lat.Ns, Np))
flush(stdout)
t = @elapsed basis = gen_basis(lat.Ns, Np)
@printf("done  (%.2f s)  RSS=%.2f GB\n", t, mem_rss_gb()); flush(stdout)

# ── 3. 只构建本节点负责的 ksector ──
println("[3/4] 构建扇区 $(seg_start)–$(seg_end) 的 ksector（保留 orbit_data + fock2rep）...")
flush(stdout)
secs_local = KSector[]
t_ks_total = @elapsed begin
    for (i, m) in enumerate(lat.ktab)
        i - 1 < seg_start && continue
        i - 1 > seg_end   && continue
        t_ks = @elapsed sec = build_ksector(basis, lat, m)
        @printf("  k=%2d  dim=%d  build %.1f s  RSS=%.2f GB\n",
                m, length(sec.reps), t_ks, mem_rss_gb())
        flush(stdout)
        push!(secs_local, sec)
    end
end
empty!(basis); GC.gc()
@printf("  ksector 完成  %d 个扇区  耗时 %.1f s  RSS=%.2f GB\n",
        length(secs_local), t_ks_total, mem_rss_gb()); flush(stdout)

# ── 预计算本节点扇区的对角元（与 φ 无关）──
println("[3b] 预计算 diag_H（φ 无关）...")
flush(stdout)
t_diag = @elapsed diag_H_local = [precompute_diag_H(s, lat, V1, V2, V3) for s in secs_local]
@printf("  完成  耗时 %.1f s  RSS=%.2f GB\n", t_diag, mem_rss_gb()); flush(stdout)

# ── 4. 遍历全部 φ 值，只算本节点的扇区 ──
println("[4/4] 谱流计算（全部 $N_phi 个 φ 点，扇区 $(seg_start)–$(seg_end)）...")
flush(stdout)

# sf_data[phi_idx] = Vector{Tuple{Int,Float64}}（k, E），只含本节点扇区的结果
sf_data = Dict{Int, Vector{Tuple{Int,Float64}}}()

for phi_idx in 0:(N_phi-1)
    phi_y = phi_grid[phi_idx + 1]   # Julia 1-based
    @printf("  φ[%2d] = %.6f (%.4f·2π) ... ", phi_idx, phi_y, phi_y/(2π))
    flush(stdout)

    t_phi = @elapsed begin
        hops = build_hops(lat, t1, t3, phi_y)

        n_local = length(secs_local)
        all_res = Vector{Vector{Tuple{Int,Float64}}}(undef, n_local)

        Threads.@threads for i in 1:n_local
            sec  = secs_local[i]
            dH   = diag_H_local[i]
            Nrep = length(sec.reps)

            if Nrep == 0
                all_res[i] = Tuple{Int,Float64}[]
                continue
            end

            nev_a = min(nev, Nrep)
            kd_a  = min(kd, Nrep)
            v0    = normalize(randn(ComplexF64, Nrep))

            function Hf(v)
                w = similar(v)
                Hv!(w, v, sec, lat, hops, V1, V2, V3, dH)
                return w
            end

            vals, _, info = eigsolve(Hf, v0, nev_a, :SR;
                ishermitian=true, krylovdim=kd_a, maxiter=500, tol=1e-10)

            info.converged < nev_a &&
                @warn "φ=$(phi_idx) k=$(sec.m) 只收敛 $(info.converged)/$nev_a 个本征值"

            n_v = min(nev_a, length(vals))
            all_res[i] = [(sec.m, e) for e in sort(real.(vals[1:n_v]))]
        end

        all_ev = Tuple{Int,Float64}[]
        for res in all_res
            append!(all_ev, res)
        end
        sort!(all_ev, by=x->x[2])
        sf_data[phi_idx] = all_ev
    end

    E0 = isempty(sf_data[phi_idx]) ? NaN : sf_data[phi_idx][1][2]
    @printf("E0=%.8f  %.1f s  RSS=%.2f GB\n", E0, t_phi, mem_rss_gb())
    flush(stdout)
end

# ── 保存结果 ──
mkpath("output_sf")
out_file = "output_sf/sf_sec$(seg_start)-$(seg_end).jld2"
jldsave(out_file;
    sf_data   = sf_data,
    phi_grid  = collect(phi_grid),
    seg_start = seg_start,
    seg_end   = seg_end,
    N_phi     = N_phi,
    Np=Np, V1=V1, V2=V2, V3=V3, t1=t1, t3=t3)
@printf("  保存: %s\n", out_file)
@printf("总耗时: %.1f s (%.2f min)\n", time()-t_total, (time()-t_total)/60)
println("完成时间: ", now())
println("完成!")
