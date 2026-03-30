# ============================================================
# Lanczos 求解器（KrylovKit.jl）
# 每个扇区流程：
#   1. precompute_diag_H  — 预计算对角势能元（一次性）
#   2. Lanczos            — 反复调用 Hv!（仅含 hopping + 预计算对角）
# 时间统计：预计算耗时 | H·v 调用次数与均耗时 | Lanczos 总耗时
# ============================================================
using KrylovKit, Printf

# ── 实际内存用量（RSS，Linux /proc，macOS fallback）──
function mem_rss_gb()
    try
        for line in eachline("/proc/self/status")
            if startswith(line, "VmRSS:")
                return parse(Int, split(line)[2]) / 1e6   # kB → GB
            end
        end
    catch
    end
    return (Sys.total_memory() - Sys.free_memory()) / 1e9  # fallback
end

function solve_sector(sec::KSector, lat::GenLat,
                      hops::Vector{Tuple{Int,Int,ComplexF64}},
                      V1::Float64, V2::Float64, V3::Float64;
                      nev::Int=4, krylovdim::Int=0,
                      verbose::Bool=false,
                      print_lock::Union{ReentrantLock,Nothing}=nothing)
    Nrep = length(sec.reps)
    Nrep == 0 && return Float64[]

    nev_actual = min(nev, Nrep)
    kd = krylovdim > 0 ? min(krylovdim, Nrep) : max(40, 2*nev_actual+20)

    # ── 预计算对角元 ──
    t_diag = @elapsed diag_H = precompute_diag_H(sec, lat, V1, V2, V3)

    v0 = normalize(randn(ComplexF64, Nrep))
    n_hv  = Ref(0)
    t_hv  = Ref(0.0)

    function Hf(v)
        result = similar(v)
        dt = @elapsed Hv!(result, v, sec, lat, hops, V1, V2, V3, diag_H)
        n_hv[] += 1
        t_hv[]  += dt
        if verbose && n_hv[] % 20 == 0
            msg = @sprintf("    [k=%2d] H·v #%3d  均 %.3f s/次  Hv累计 %.1f s\n",
                           sec.m, n_hv[], t_hv[]/n_hv[], t_hv[])
            if print_lock !== nothing
                lock(print_lock) do; print(msg); flush(stdout); end
            else
                print(msg); flush(stdout)
            end
        end
        return result
    end

    t_eig = @elapsed vals, _, info = eigsolve(Hf, v0, nev_actual, :SR;
        ishermitian=true, krylovdim=kd, maxiter=500, tol=1e-10)

    info.converged < nev_actual &&
        @warn "k=$(sec.m) 只收敛 $(info.converged)/$nev_actual 个本征值"

    # ── 本扇区耗时明细 ──
    nhv = n_hv[]
    msg = @sprintf(
        "    [k=%2d] 预计算对角 %5.2f s | H·v ×%3d 共 %6.1f s (均 %.3f s/次) | Lanczos %6.1f s\n",
        sec.m, t_diag, nhv, t_hv[], nhv > 0 ? t_hv[]/nhv : 0.0, t_eig)
    if print_lock !== nothing
        lock(print_lock) do; print(msg); flush(stdout); end
    else
        print(msg); flush(stdout)
    end

    return sort(real.(vals[1:min(nev_actual, length(vals))]))
end

function solve_sector_with_vec(sec::KSector, lat::GenLat,
                               hops::Vector{Tuple{Int,Int,ComplexF64}},
                               V1::Float64, V2::Float64, V3::Float64;
                               nev::Int=1, krylovdim::Int=0,
                               verbose::Bool=false,
                               print_lock::Union{ReentrantLock,Nothing}=nothing)
    Nrep = length(sec.reps)
    Nrep == 0 && return Float64[], Vector{ComplexF64}[]

    nev_actual = min(nev, Nrep)
    kd = krylovdim > 0 ? min(krylovdim, Nrep) : max(40, 2*nev_actual+20)

    # ── 预计算对角元 ──
    t_diag = @elapsed diag_H = precompute_diag_H(sec, lat, V1, V2, V3)

    v0 = normalize(randn(ComplexF64, Nrep))
    n_hv = Ref(0)
    t_hv = Ref(0.0)

    function Hf(v)
        result = similar(v)
        dt = @elapsed Hv!(result, v, sec, lat, hops, V1, V2, V3, diag_H)
        n_hv[] += 1
        t_hv[]  += dt
        if verbose && n_hv[] % 20 == 0
            msg = @sprintf("    [k=%2d] H·v #%3d  均 %.3f s/次  Hv累计 %.1f s\n",
                           sec.m, n_hv[], t_hv[]/n_hv[], t_hv[])
            if print_lock !== nothing
                lock(print_lock) do; print(msg); flush(stdout); end
            else
                print(msg); flush(stdout)
            end
        end
        return result
    end

    t_eig = @elapsed vals, vecs, info = eigsolve(Hf, v0, nev_actual, :SR;
        ishermitian=true, krylovdim=kd, maxiter=500, tol=1e-10)

    info.converged < nev_actual &&
        @warn "k=$(sec.m) 只收敛 $(info.converged)/$nev_actual 个本征值"

    nhv = n_hv[]
    msg = @sprintf(
        "    [k=%2d] 预计算对角 %5.2f s | H·v ×%3d 共 %6.1f s (均 %.3f s/次) | Lanczos %6.1f s\n",
        sec.m, t_diag, nhv, t_hv[], nhv > 0 ? t_hv[]/nhv : 0.0, t_eig)
    if print_lock !== nothing
        lock(print_lock) do; print(msg); flush(stdout); end
    else
        print(msg); flush(stdout)
    end

    n = min(nev_actual, length(vals))
    return real.(vals[1:n]), vecs[1:n]
end

function compute_spectrum(secs::Vector{KSector}, lat::GenLat,
                          hops::Vector{Tuple{Int,Int,ComplexF64}},
                          V1::Float64, V2::Float64, V3::Float64;
                          nev::Int=4, krylovdim::Int=0, verbose::Bool=false)
    n  = length(secs)
    all_res = Vector{Vector{Tuple{Int,Float64}}}(undef, n)
    lk = ReentrantLock()

    Threads.@threads for i in 1:n
        sec   = secs[i]
        t_sec = @elapsed evals = solve_sector(sec, lat, hops, V1, V2, V3;
                                               nev=nev, krylovdim=krylovdim,
                                               verbose=verbose, print_lock=lk)
        all_res[i] = [(sec.m, e) for e in evals]
        lock(lk) do
            isempty(evals) || @printf(
                "  ✓ k=%2d  dim=%7d  E_min=%10.6f  扇区总耗时 %6.1f s\n",
                sec.m, length(sec.reps), evals[1], t_sec)
            flush(stdout)
        end
    end

    all_ev = Tuple{Int,Float64}[]
    for res in all_res
        append!(all_ev, res)
    end
    sort!(all_ev, by=x->x[2])
    return all_ev
end

# 同 compute_spectrum，但同时保留每个扇区的基态波函数。
# 返回 (all_ev, gs_vecs)，其中 gs_vecs::Dict{Int,Vector{ComplexF64}} 是 k → 基态向量。
# 避免后处理阶段重复 Lanczos。
function compute_spectrum_with_vecs(secs::Vector{KSector}, lat::GenLat,
                                    hops::Vector{Tuple{Int,Int,ComplexF64}},
                                    V1::Float64, V2::Float64, V3::Float64;
                                    nev::Int=4, krylovdim::Int=0, verbose::Bool=false)
    n        = length(secs)
    all_res  = Vector{Vector{Tuple{Int,Float64}}}(undef, n)
    all_vecs = Vector{Vector{ComplexF64}}(undef, n)
    lk = ReentrantLock()

    Threads.@threads for i in 1:n
        sec = secs[i]
        t_sec = @elapsed begin
            evals, vecs = solve_sector_with_vec(sec, lat, hops, V1, V2, V3;
                                                nev=nev, krylovdim=krylovdim,
                                                verbose=verbose, print_lock=lk)
        end
        all_res[i]  = [(sec.m, e) for e in evals]
        all_vecs[i] = isempty(vecs) ? ComplexF64[] : vecs[1]   # 取最低本征态
        lock(lk) do
            isempty(evals) || @printf(
                "  ✓ k=%2d  dim=%7d  E_min=%10.6f  扇区总耗时 %6.1f s\n",
                sec.m, length(sec.reps), evals[1], t_sec)
            flush(stdout)
        end
    end

    all_ev = Tuple{Int,Float64}[]
    for res in all_res
        append!(all_ev, res)
    end
    sort!(all_ev, by=x->x[2])

    gs_vecs = Dict{Int,Vector{ComplexF64}}(
        secs[i].m => all_vecs[i] for i in 1:n if !isempty(all_vecs[i])
    )
    return all_ev, gs_vecs
end

# ============================================================
# CSR 稀疏矩阵求解器
# 流程：
#   A. Threads.@threads 并行 build_sparse_H（一次性，fock2rep 查找在此完成）
#   B. Threads.@threads 并行 Lanczos（mul!(w,H,v)，顺序内存访问）
#   C. 返回 (all_ev, gs_vecs)，供后处理直接使用
# ============================================================

"""
    compute_spectrum_sparse_with_vecs(secs, lat, hops, V1, V2, V3; ...)

【旧接口，保留供兼容】批量建 CSR + 并行 Lanczos。
内存峰值高（所有扇区的 orbit_data/fock2rep 同时存在）。
新代码请在 main.jl 中用流水线方式手动调用 build_sparse_H + lanczos_sparse_sectors。
"""
function compute_spectrum_sparse_with_vecs(
        secs::Vector{KSector}, lat::GenLat,
        hops::Vector{Tuple{Int,Int,ComplexF64}},
        V1::Float64, V2::Float64, V3::Float64;
        nev::Int=4, krylovdim::Int=120, verbose::Bool=false)

    n  = length(secs)
    lk = ReentrantLock()

    # ── A. 并行建 CSR ──
    @printf("  [CSR-A] 并行构建 %d 个稀疏 H 矩阵 (%d 线程)...\n",
            n, Threads.nthreads())
    @printf("  [MEM] CSR 构建前 RSS = %.2f GB\n", mem_rss_gb())
    flush(stdout)
    H_csrs = Vector{SparseMatrixCSC{ComplexF64,Int32}}(undef, n)
    t_build = @elapsed Threads.@threads for i in 1:n
        H_csrs[i] = build_sparse_H(secs[i], lat, hops, V1, V2, V3)
        empty!(secs[i].fock2rep)
        lock(lk) do
            @printf("  [MEM] k=%2d CSR 完成 nnz=%d (%.2f GB)  RSS=%.2f GB\n",
                    secs[i].m, nnz(H_csrs[i]),
                    nnz(H_csrs[i])*12/1e9, mem_rss_gb())
            flush(stdout)
        end
    end
    nnz_total = sum(nnz(H) for H in H_csrs)
    @printf("  [CSR-A] 完成  耗时 %.1f s  总非零元 %d  内存估算 %.2f GB\n",
            t_build, nnz_total, nnz_total * 12 / 1e9)
    @printf("  [MEM] CSR 全部完成后 RSS = %.2f GB\n", mem_rss_gb())
    flush(stdout)

    return lanczos_sparse_sectors(secs, H_csrs; nev=nev, krylovdim=krylovdim, verbose=verbose)
end

"""
    lanczos_sparse_sectors(secs, H_csrs; nev, krylovdim, verbose)

对预建的 CSR 矩阵并行 Lanczos 求解。
secs 只需保留 .m 和 .reps（orbit_data/fock2rep 可已清空）。
供流水线模式的 main.jl 调用，避免 orbit_data/fock2rep 在 Lanczos 阶段仍占用内存。
"""
function lanczos_sparse_sectors(
        secs::Vector{KSector},
        H_csrs::Vector{SparseMatrixCSC{ComplexF64,Int32}};
        nev::Int=4, krylovdim::Int=60, verbose::Bool=false)

    n  = length(secs)
    lk = ReentrantLock()

    @printf("  [CSR-B] 并行 Lanczos  krylovdim=%d  nev=%d  %d 个扇区...\n",
            krylovdim, nev, n)
    @printf("  [MEM] Lanczos 开始前 RSS = %.2f GB  (预估 Krylov 峰值 +%.1f GB)\n",
            mem_rss_gb(), n * krylovdim * (length(secs[1].reps)*16/1e9))
    flush(stdout)
    all_res  = Vector{Vector{Tuple{Int,Float64}}}(undef, n)
    all_vecs = Vector{Vector{ComplexF64}}(undef, n)

    t_solve = @elapsed Threads.@threads for i in 1:n
        sec  = secs[i]
        H    = H_csrs[i]
        Nrep = length(sec.reps)

        if Nrep == 0
            all_res[i]  = Tuple{Int,Float64}[]
            all_vecs[i] = ComplexF64[]
            continue
        end

        nev_actual = min(nev, Nrep)
        kd  = min(krylovdim, Nrep)
        v0  = normalize(randn(ComplexF64, Nrep))
        n_hv = Ref(0); t_hv = Ref(0.0)

        function Hf(v)
            result = similar(v)
            dt = @elapsed mul!(result, H, v)
            n_hv[] += 1; t_hv[] += dt
            if verbose && n_hv[] % 50 == 0
                msg = @sprintf("    [k=%2d] H·v #%3d  均 %.5f s/次\n",
                               sec.m, n_hv[], t_hv[] / n_hv[])
                lock(lk) do; print(msg); flush(stdout); end
            end
            return result
        end

        t_eig = @elapsed vals, vecs, info = eigsolve(
            Hf, v0, nev_actual, :SR;
            ishermitian=true, krylovdim=kd, maxiter=500, tol=1e-10)

        info.converged < nev_actual &&
            @warn "k=$(sec.m) 只收敛 $(info.converged)/$nev_actual 个本征值"

        nhv = n_hv[]
        msg = @sprintf(
            "  ✓ [k=%2d] dim=%7d  H·v×%3d  均 %.5f s/次  Lanczos %.1f s  E0=%.8f  RSS=%.2f GB\n",
            sec.m, Nrep, nhv, nhv > 0 ? t_hv[] / nhv : 0.0, t_eig,
            isempty(vals) ? NaN : real(vals[1]), mem_rss_gb())
        lock(lk) do; print(msg); flush(stdout); end

        n_vals = min(nev_actual, length(vals))
        all_res[i]  = [(sec.m, e) for e in sort(real.(vals[1:n_vals]))]
        all_vecs[i] = isempty(vecs) ? ComplexF64[] : vecs[1]
    end

    @printf("  [CSR-B] Lanczos 完成  耗时 %.1f s (%.2f min)\n",
            t_solve, t_solve / 60)
    @printf("  [MEM] Lanczos 完成后 RSS = %.2f GB\n", mem_rss_gb())
    flush(stdout)

    all_ev = Tuple{Int,Float64}[]
    for res in all_res
        append!(all_ev, res)
    end
    sort!(all_ev, by=x->x[2])

    gs_vecs = Dict{Int,Vector{ComplexF64}}(
        secs[i].m => all_vecs[i]
        for i in 1:n if isassigned(all_vecs, i) && !isempty(all_vecs[i])
    )
    return all_ev, gs_vecs
end
