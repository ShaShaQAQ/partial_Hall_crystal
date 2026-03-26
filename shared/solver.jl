# ============================================================
# Lanczos 求解器（KrylovKit.jl）
# 每个扇区流程：
#   1. precompute_diag_H  — 预计算对角势能元（一次性）
#   2. Lanczos            — 反复调用 Hv!（仅含 hopping + 预计算对角）
# 时间统计：预计算耗时 | H·v 调用次数与均耗时 | Lanczos 总耗时
# ============================================================
using KrylovKit, Printf

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
