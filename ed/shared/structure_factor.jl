# ============================================================
# 密度结构因子 N(k)（通用版本）
#
# 采用原胞归一化：
#   N(k) = (1/Nuc) <0| Σ_{j,l} e^{ik(Rj-Rl)} (n̂j-n̄)(n̂l-n̄) |0>
#
# 等价写法：N(k) = (1/Nuc) <|ρ_k^-|²>
#   ρ_k^-(F) = Σ_{i∈occ(F)} e^{ik·ri}  -  n̄ · Σ_i e^{ik·ri}
#
# 注：对团簇 k≠0，Σ_i e^{ik·ri} = 0（正交性），结果同旧公式。
#     k=0 时 ρ_0^- = Np - n̄*Ns = 0，故 n(0) = 0。
# ============================================================

using Printf

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

"""
    structure_factor_survey_representatives(sec, v, lat, Np)

Compute the unit-cell-normalized structure factor at all allowed torus
momenta without expanding the momentum state into the full Fock basis.
For these momenta, `abs2(rho_q)` is invariant under every translation in a
representative's orbit, so a diagonal translation-invariant observable can
be accumulated once per representative with weight `abs2(v[i])`.
"""
function structure_factor_survey_representatives(
        sec::KSector,
        v::Vector{ComplexF64},
        lat::GenLat,
        Np::Int)
    length(v) == length(sec.reps) ||
        throw(DimensionMismatch("state does not match momentum sector"))
    momentum_count = length(lat.kpoints)
    phases = Matrix{ComplexF64}(undef, momentum_count, lat.Ns)
    phase_sums = zeros(ComplexF64, momentum_count)
    for (q_index, momentum) in enumerate(lat.kpoints)
        for (site, (ix, iy)) in enumerate(lat.sites)
            rx, ry = site_uc_pos(lat, ix, iy)
            phase = cis(momentum[1] * rx + momentum[2] * ry)
            phases[q_index, site] = phase
            phase_sums[q_index] += phase
        end
    end

    thread_totals = [zeros(Float64, momentum_count)
                     for _ in 1:Threads.nthreads()]
    rho_buffers = [zeros(ComplexF64, momentum_count)
                   for _ in 1:Threads.nthreads()]
    filling = Np / lat.Ns
    Threads.@threads :static for representative_index in eachindex(sec.reps)
        tid = Threads.threadid()
        rho = rho_buffers[tid]
        fill!(rho, 0.0 + 0.0im)
        fock_state = sec.reps[representative_index]
        occupied = fock_state
        while occupied != 0
            site = trailing_zeros(occupied) + 1
            @inbounds for q_index in 1:momentum_count
                rho[q_index] += phases[q_index, site]
            end
            occupied &= occupied - 1
        end

        weight = abs2(v[representative_index])
        totals = thread_totals[tid]
        @inbounds for q_index in 1:momentum_count
            centered_density = rho[q_index] -
                               filling * phase_sums[q_index]
            totals[q_index] += weight * abs2(centered_density)
        end
    end

    norm2 = sum(abs2, v)
    norm2 > 0 || throw(ArgumentError("state has zero norm"))
    structure_factors = zeros(Float64, momentum_count)
    for totals in thread_totals
        structure_factors .+= totals
    end
    structure_factors ./= lat.Nuc * norm2
    particle_number = sum(
        abs2(v[i]) * count_ones(sec.reps[i]) for i in eachindex(sec.reps)
    ) / norm2

    results_ordered = Tuple{Int,Float64,Float64,Float64}[]
    for (index, (m, momentum)) in enumerate(zip(lat.ktab, lat.kpoints))
        push!(results_ordered,
              (m, structure_factors[index], momentum[1], momentum[2]))
    end
    results_sorted = sort(results_ordered, by=result -> -result[2])
    return results_ordered, results_sorted, norm2, particle_number
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
