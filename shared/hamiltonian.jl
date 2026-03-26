# ============================================================
# Matrix-free 哈密顿量：H|v> 在动量扇区中的作用
#
# 优化：precompute_diag_H 预先算好每个代表元的对角势能元，
#       Hv! 中直接用，省去每次 Lanczos 步的 Ns×6 嵌套循环
#       和 fock2rep 哈希查找（对角项不需要查表）。
# ============================================================

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
    for i in 1:Nrep
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
             diag_H::Vector{Float64})
    Nrep = length(sec.reps)
    fill!(w, 0.0+0im)

    for i in 1:Nrep
        abs2(v[i]) < 1e-20 && continue
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

    Is = Int32[]
    Js = Int32[]
    Vs = ComplexF64[]
    sizehint!(Is, Nrep * 50)
    sizehint!(Js, Nrep * 50)
    sizehint!(Vs, Nrep * 50)

    for i in 1:Nrep
        ni  = 1.0 / sec.norms[i]
        ni2 = ni * ni

        # ── 对角势能（内联 precompute_diag_H） ──
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
        push!(Is, Int32(i)); push!(Js, Int32(i)); push!(Vs, ComplexF64(diag_val))

        # ── 跳跃（off-diagonal） ──
        for (c_α, F) in sec.orbit_data[i]
            coeff_base = c_α * ni
            for (tgt, src, tval) in hops
                (F>>(src-1))&1 == 0 && continue
                (F>>(tgt-1))&1 == 1 && continue
                newF, sgn_hop = hop(F, tgt, src)
                info = get(sec.fock2rep, newF, nothing)
                info === nothing && continue
                (j, c_conj_norm) = info
                push!(Is, Int32(j))
                push!(Js, Int32(i))
                push!(Vs, coeff_base * tval * sgn_hop * c_conj_norm)
            end
        end
    end

    # Julia sparse() 自动合并重复 (i,j) 索引（求和），正确处理多轨道贡献
    return sparse(Is, Js, Vs, Nrep, Nrep)
end
