# ============================================================
# 密度结构因子 N(k)（通用版本）
#
# 论文 Kourtis et al. PRB 86, 235118 (2012) Eq.(27)：
#   n(k) = (1/Ns) <0| Σ_{j,l} e^{ik(Rj-Rl)} (n̂j-n̄)(n̂l-n̄) |0>
#
# 等价写法：n(k) = (1/Ns) <|ρ_k^-|²>
#   ρ_k^-(F) = Σ_{i∈occ(F)} e^{ik·ri}  -  n̄ · Σ_i e^{ik·ri}
#
# 注：对团簇 k≠0，Σ_i e^{ik·ri} = 0（正交性），结果同旧公式。
#     k=0 时 ρ_0^- = Np - n̄*Ns = 0，故 n(0) = 0。
# ============================================================

# 格点物理坐标：(ix,iy) → ix*a1 + iy*a2。
function site_phys_pos(lat::GenLat, ix::Int, iy::Int)
    rx = Float64(ix) * lat.a1[1] + Float64(iy) * lat.a2[1]
    ry = Float64(ix) * lat.a1[2] + Float64(iy) * lat.a2[2]
    return rx, ry
end

# 原胞 Bravais 格矢：同一原胞的两个子格位点返回相同坐标
# 用于 Eq.(27) 的 Rj，保证 n(q) 对完整倒格矢 b1,b2 周期
function site_uc_pos(lat::GenLat, ix::Int, iy::Int)
    uc_y = (iy - mod(iy, 2)) ÷ 2
    rx = Float64(ix) * lat.a1[1] + Float64(uc_y) * 2.0 * lat.a2[1]
    ry = Float64(ix) * lat.a1[2] + Float64(uc_y) * 2.0 * lat.a2[2]
    return rx, ry
end

# 将 K-扇区本征态展开到 Fock 空间
function expand_to_fock(sec::KSector, v::Vector{ComplexF64})
    coeffs = Dict{Int64,ComplexF64}()
    sizehint!(coeffs, length(sec.fock2rep))
    for (F, (rep_idx, c_conj_norm)) in sec.fock2rep
        coeffs[F] = v[rep_idx] * conj(c_conj_norm)
    end
    return coeffs
end

# 在给定波矢 (qx,qy) 处计算 n(k)，严格按 Eq.(27) 减去 n̄ 项
# n̄ = Np / Ns
function compute_sq(coeffs::Dict{Int64,ComplexF64},
                    lat::GenLat, qx::Float64, qy::Float64, n̄::Float64)
    Ns = lat.Ns
    eiqr = Vector{ComplexF64}(undef, Ns)
    phase_sum = zero(ComplexF64)          # Σ_i e^{iq·ri}
    for (s, (ix, iy)) in enumerate(lat.sites)
        rx, ry = site_uc_pos(lat, ix, iy) # 用原胞 Bravais 位置，保证 n(q+G)=n(q)
        eiqr[s]   = cis(qx*rx + qy*ry)
        phase_sum += eiqr[s]
    end

    norm2 = sum(abs2(c) for c in values(coeffs))
    Sq = 0.0
    for (F, c_F) in coeffs
        rho_q = zero(ComplexF64)
        tmp = F
        while tmp != 0
            s = trailing_zeros(tmp) + 1
            rho_q += eiqr[s]
            tmp &= tmp - 1
        end
        rho_q_minus = rho_q - n̄ * phase_sum   # Eq.(27): subtract n̄ term
        Sq += abs2(c_F) * abs2(rho_q_minus)
    end
    return Sq / (lat.Nuc * norm2)
end

# ── 主函数：计算所有 Nuc 个团簇 k 点处的 n(k) ──
function structure_factor_survey(sec::KSector, v::Vector{ComplexF64},
                                 lat::GenLat, Np::Int)
    coeffs = expand_to_fock(sec, v)

    norm2    = sum(abs2(c) for c in values(coeffs))
    Np_check = sum(abs2(c)*count_ones(F) for (F,c) in coeffs) / norm2

    n̄ = Np / lat.Ns

    results_ordered = Tuple{Int,Float64,Float64,Float64}[]
    for (m, kphys) in zip(lat.ktab, lat.kpoints)
        sq = compute_sq(coeffs, lat, kphys[1], kphys[2], n̄)
        push!(results_ordered, (m, sq, kphys[1], kphys[2]))
    end
    results_sorted = sort(results_ordered, by=x->-x[2])

    return results_ordered, results_sorted, norm2, Np_check
end

# 打印报告
function print_sq_report(results_sorted, norm2, Np_check, Np, Ns, k_sector)
    println("\n" * "="^60)
    println("结构因子 n(k) [Eq.27] — 基态扇区 k=$(k_sector)")
    @printf("  归一化检验: Σ|c_F|² = %.6f\n", norm2)
    @printf("  粒子数检验: <N> = %.4f (应为 %d)\n", Np_check, Np)
    println("-"^60)
    println("  排名   ik    n(k)      (kx, ky)")
    println("-"^60)
    for (i, (m, sq, kx, ky)) in enumerate(results_sorted[1:min(10,end)])
        @printf("  %3d.   %2d   %8.4f   (%.4f, %.4f)\n", i, m, sq, kx, ky)
    end
    println("="^60)
end
