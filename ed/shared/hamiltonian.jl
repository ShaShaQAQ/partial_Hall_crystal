# ============================================================
# Matrix-free 哈密顿量：H|v> 在动量扇区中的作用
#
# 优化：precompute_diag_H 预先算好每个代表元的对角势能元，
#       Hv! 中直接用，省去每次 Lanczos 步的 Ns×6 嵌套循环
#       和 fock2rep 哈希查找（对角项不需要查表）。
# ============================================================
using Printf
# mem_rss_gb 在 solver.jl 中定义；若单独加载本文件则用此备用定义
if !isdefined(Main, :mem_rss_gb)
    function mem_rss_gb()
        for line in eachline("/proc/self/status")
            if startswith(line, "VmRSS:")
                return parse(Int, split(line)[2]) / 1e6
            end
        end
        return NaN
    end
end

"""
    precompute_diag_H(sec, lat, V1, V2, V3) -> Vector{Float64}

预计算动量扇区内每个代表元 i 的对角哈密顿量元（势能部分）：

    diag_H[i] = (1/||i||²) Σ_{(c_α,F)∈orbit_i} |c_α|² · Eint(F)

调用一次，之后每次 Hv! 直接 w[i] += v[i]*diag_H[i]，
省去 O(Nrep × orbit_size × Ns × 6) 的重复计算。
"""
function precompute_diag_H(sec::KSector, lat::GenLat,
                            V1::Float64, V2::Float64, V3::Float64)
    Ns   = lat.Ns
    Nrep = length(sec.reps)
    diag_H = zeros(Float64, Nrep)
    Threads.@threads :static for i in 1:Nrep
        ni2 = 1.0 / sec.norms[i]^2
        for (c_α, F) in sec.orbit_data[i]
            Eint = 0.0
            for s in 1:Ns
                (F>>(s-1))&1 == 0 && continue
                for k in 1:6
                    Eint += V1 * ((F>>(lat.nb1[s,k]-1))&1) * 0.5
                    Eint += V2 * ((F>>(lat.nb2[s,k]-1))&1) * 0.5
                    Eint += V3 * ((F>>(lat.nb3[s,k]-1))&1) * 0.5
                end
            end
            diag_H[i] += abs2(c_α) * ni2 * Eint
        end
    end
    return diag_H
end

# ── 带预计算对角元的版本（Lanczos 求解时使用）──
function Hv!(w::AbstractVector{ComplexF64}, v::AbstractVector{ComplexF64},
             sec::KSector, lat::GenLat,
             hops::Vector{Tuple{Int,Int,ComplexF64}},
             V1::Float64, V2::Float64, V3::Float64,
             diag_H::Vector{Float64};
             skip_tol::Float64=1e-20)
    Nrep = length(sec.reps)
    fill!(w, 0.0+0im)

    for i in 1:Nrep
        abs2(v[i]) < skip_tol && continue
        vi = v[i]
        ni = 1.0 / sec.norms[i]

        # 对角势能（O(1)，预计算值）
        @inbounds w[i] += vi * diag_H[i]

        for (c_α, F) in sec.orbit_data[i]
            coeff_F = vi * c_α * ni

            # 动能：hopping
            for (tgt, src, tval) in hops
                (F>>(src-1))&1 == 0 && continue
                (F>>(tgt-1))&1 == 1 && continue
                newF, sgn_hop = hop(F, tgt, src)
                info = get(sec.fock2rep, newF, nothing)
                info === nothing && continue
                (j, c_conj_norm) = info
                @inbounds w[j] += coeff_F * tval * sgn_hop * c_conj_norm
            end
        end
    end
    return w
end

"""Allocate one scatter output buffer per Julia thread for `Hv_threaded!`."""
function threaded_Hv_workspace(Nrep::Int)
    return [zeros(ComplexF64, Nrep) for _ in 1:Threads.nthreads()]
end

function Hv_threaded!(
        w::AbstractVector{ComplexF64},
        v::AbstractVector{ComplexF64},
        sec::KSector,
        lat::GenLat,
        hops::Vector{Tuple{Int,Int,ComplexF64}},
        diag_H::Vector{Float64},
        scratch::Vector{Vector{ComplexF64}};
        skip_tol::Float64=1e-20)
    Nrep = length(sec.reps)
    length(w) == Nrep || throw(DimensionMismatch("output has wrong length"))
    length(v) == Nrep || throw(DimensionMismatch("input has wrong length"))
    length(diag_H) == Nrep ||
        throw(DimensionMismatch("diagonal has wrong length"))
    length(scratch) == Threads.nthreads() ||
        throw(DimensionMismatch(
            "workspace must have one vector per Julia thread"))
    all(length(buffer) == Nrep for buffer in scratch) ||
        throw(DimensionMismatch("workspace vector has wrong length"))

    Threads.@threads :static for tid in eachindex(scratch)
        fill!(scratch[tid], 0.0 + 0.0im)
    end

    Threads.@threads :static for i in 1:Nrep
        abs2(v[i]) < skip_tol && continue
        wi = scratch[Threads.threadid()]
        vi = v[i]
        ni = 1.0 / sec.norms[i]
        @inbounds wi[i] += vi * diag_H[i]

        for (coefficient, state) in sec.orbit_data[i]
            state_coefficient = vi * coefficient * ni
            for (target, source, hopping) in hops
                (state >> (source - 1)) & 1 == 0 && continue
                (state >> (target - 1)) & 1 == 1 && continue
                new_state, hopping_sign = hop(state, target, source)
                info = get(sec.fock2rep, new_state, nothing)
                info === nothing && continue
                j, conjugate_coefficient = info
                @inbounds wi[j] += state_coefficient * hopping *
                                           hopping_sign * conjugate_coefficient
            end
        end
    end

    Threads.@threads :static for i in 1:Nrep
        value = 0.0 + 0.0im
        @inbounds for tid in eachindex(scratch)
            value += scratch[tid][i]
        end
        w[i] = value
    end
    return w
end

# ── 原始版本（不依赖预计算，兼容旧调用）──
function Hv!(w::AbstractVector{ComplexF64}, v::AbstractVector{ComplexF64},
             sec::KSector, lat::GenLat,
             hops::Vector{Tuple{Int,Int,ComplexF64}},
             V1::Float64, V2::Float64, V3::Float64)
    Ns = lat.Ns
    Nrep = length(sec.reps)
    fill!(w, 0.0+0im)

    for i in 1:Nrep
        abs2(v[i]) < 1e-20 && continue
        vi = v[i]
        ni = 1.0 / sec.norms[i]

        for (c_α, F) in sec.orbit_data[i]
            coeff_F = vi * c_α * ni

            # 动能：hopping
            for (tgt, src, tval) in hops
                (F>>(src-1))&1 == 0 && continue
                (F>>(tgt-1))&1 == 1 && continue
                newF, sgn_hop = hop(F, tgt, src)
                info = get(sec.fock2rep, newF, nothing)
                info === nothing && continue
                (j, c_conj_norm) = info
                @inbounds w[j] += coeff_F * tval * sgn_hop * c_conj_norm
            end

            # 势能：对角项（实时计算）
            Eint = 0.0
            for s in 1:Ns
                (F>>(s-1))&1 == 0 && continue
                for k in 1:6
                    Eint += V1 * ((F>>(lat.nb1[s,k]-1))&1) * 0.5
                    Eint += V2 * ((F>>(lat.nb2[s,k]-1))&1) * 0.5
                    Eint += V3 * ((F>>(lat.nb3[s,k]-1))&1) * 0.5
                end
            end
            if Eint != 0.0
                info = get(sec.fock2rep, F, nothing)
                if info !== nothing
                    (j, c_conj_norm) = info
                    @inbounds w[j] += coeff_F * Eint * c_conj_norm
                end
            end
        end
    end
    return w
end

# ============================================================
# CSR 预建稀疏哈密顿量
#
# 一次性将 Hv! 的所有操作展开为 SparseMatrixCSC 的非零元，
# 后续 Hv = mul!(w, H, v)，顺序内存访问，BLAS 优化。
# ============================================================
using SparseArrays

"""
    build_sparse_H(sec, lat, hops, V1, V2, V3) -> SparseMatrixCSC{ComplexF64,Int32}

将动量扇区 sec 的哈密顿量预建为稀疏矩阵（CSR 格式）。
包含：对角势能（内联 precompute_diag_H 逻辑）+ off-diagonal 跳跃项。
fock2rep 哈希查找只发生一次（此函数内），之后 Lanczos 全程无 Dict 访问。
"""
function build_sparse_H(sec::KSector, lat::GenLat,
                        hops::Vector{Tuple{Int,Int,ComplexF64}},
                        V1::Float64, V2::Float64, V3::Float64)
    Nrep = length(sec.reps)
    Nrep == 0 && return spzeros(ComplexF64, Int32, 0, 0)
    Ns = lat.Ns

    @printf("    [build_sparse_H k=%2d] 开始 Nrep=%d  RSS=%.2f GB\n",
            sec.m, Nrep, mem_rss_gb()); flush(stdout)

    # ── 第一遍：统计每列唯一非零行数，用于精确分配 CSC 数组 ──
    col_nnz = zeros(Int32, Nrep)
    row_sets = [Set{Int32}() for _ in 1:Threads.nthreads()]
    for rows in row_sets
        sizehint!(rows, 256)
    end
    Threads.@threads :static for i in 1:Nrep
        row_set = row_sets[Threads.threadid()]
        empty!(row_set)
        push!(row_set, Int32(i))                    # 对角元
        for (_, F) in sec.orbit_data[i]
            for (tgt, src, _) in hops
                (F>>(src-1))&1 == 0 && continue
                (F>>(tgt-1))&1 == 1 && continue
                newF, _ = hop(F, tgt, src)
                info = get(sec.fock2rep, newF, nothing)
                info === nothing && continue
                push!(row_set, info[1])
            end
        end
        col_nnz[i] = Int32(length(row_set))
    end

    colptr = Vector{Int32}(undef, Nrep + 1)
    colptr[1] = Int32(1)
    for i in 1:Nrep
        colptr[i+1] = colptr[i] + col_nnz[i]
    end
    total_nnz = Int(colptr[Nrep+1]) - 1

    @printf("    [build_sparse_H k=%2d] 第一遍完成  nnz=%d  RSS=%.2f GB\n",
            sec.m, total_nnz, mem_rss_gb()); flush(stdout)

    rowval = Vector{Int32}(undef, total_nnz)
    nzval  = Vector{ComplexF64}(undef, total_nnz)

    @printf("    [build_sparse_H k=%2d] CSC 数组分配后  RSS=%.2f GB\n",
            sec.m, mem_rss_gb()); flush(stdout)

    # ── 第二遍：逐列用小 Dict 累加，直接写入 CSC 数组 ──
    col_buffers = [Dict{Int32,ComplexF64}()
                   for _ in 1:Threads.nthreads()]
    column_keys = [Int32[] for _ in 1:Threads.nthreads()]
    for buffer in col_buffers
        sizehint!(buffer, 256)
    end
    for keys_buffer in column_keys
        sizehint!(keys_buffer, 256)
    end

    Threads.@threads :static for i in 1:Nrep
        tid = Threads.threadid()
        col_buf = col_buffers[tid]
        keys_buffer = column_keys[tid]
        ni  = 1.0 / sec.norms[i]
        ni2 = ni * ni
        empty!(col_buf)

        # 对角势能
        diag_val = 0.0
        for (c_α, F) in sec.orbit_data[i]
            Eint = 0.0
            for s in 1:Ns
                (F>>(s-1))&1 == 0 && continue
                for k in 1:6
                    Eint += V1 * ((F>>(lat.nb1[s,k]-1))&1) * 0.5
                    Eint += V2 * ((F>>(lat.nb2[s,k]-1))&1) * 0.5
                    Eint += V3 * ((F>>(lat.nb3[s,k]-1))&1) * 0.5
                end
            end
            diag_val += abs2(c_α) * ni2 * Eint
        end
        col_buf[Int32(i)] = ComplexF64(diag_val)

        # 跳跃（off-diagonal）
        for (c_α, F) in sec.orbit_data[i]
            coeff_base = c_α * ni
            for (tgt, src, tval) in hops
                (F>>(src-1))&1 == 0 && continue
                (F>>(tgt-1))&1 == 1 && continue
                newF, sgn_hop = hop(F, tgt, src)
                info = get(sec.fock2rep, newF, nothing)
                info === nothing && continue
                (j, c_conj_norm) = info
                v = coeff_base * tval * sgn_hop * c_conj_norm
                col_buf[j] = get(col_buf, j, zero(ComplexF64)) + v
            end
        end

        # 写入 CSC（行号排序，Julia SpMV 要求）
        empty!(keys_buffer)
        append!(keys_buffer, keys(col_buf))
        sort!(keys_buffer)
        pos = Int(colptr[i])
        for j in keys_buffer
            rowval[pos] = j
            nzval[pos]  = col_buf[j]
            pos += 1
        end
    end

    @printf("    [build_sparse_H k=%2d] 第二遍完成  RSS=%.2f GB\n",
            sec.m, mem_rss_gb()); flush(stdout)

    return SparseMatrixCSC(Nrep, Nrep, colptr, rowval, nzval)
end
