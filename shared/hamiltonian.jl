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
