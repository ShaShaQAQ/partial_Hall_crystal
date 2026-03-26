# ============================================================
# 实空间波函数分析模块
#
# 依赖：structure_factor.jl 中的 site_phys_pos、expand_to_fock
#
# 提供：
#   site_occupation      — 格点占据数 ⟨n_i⟩
#   structure_factor_at  — 单个 q 点处的 n(q)，按 Kourtis Eq.(27)
#   structure_factor_bz  — 均匀 BZ 网格扫描 n(q)
# ============================================================

# ── 格点占据数 ⟨n_i⟩ ──
# 返回长度为 Ns 的向量，∑_i ⟨n_i⟩ = Np
function site_occupation(coeffs::Dict{Int64,ComplexF64}, lat::GenLat)
    occ = zeros(Float64, lat.Ns)
    norm2 = sum(abs2(c) for c in values(coeffs))
    for (F, c) in coeffs
        w = abs2(c)
        tmp = F
        while tmp != 0
            s = trailing_zeros(tmp) + 1
            occ[s] += w
            tmp &= tmp - 1
        end
    end
    occ ./= norm2
    return occ
end

# ── 单个 q 点处的 n(q)，严格按 Kourtis et al. Eq.(27) ──
# n̄ = Np / Ns
function structure_factor_at(coeffs::Dict{Int64,ComplexF64},
                              lat::GenLat, qx::Float64, qy::Float64, n̄::Float64)
    Ns = lat.Ns
    eiqr = Vector{ComplexF64}(undef, Ns)
    phase_sum = zero(ComplexF64)           # Σ_i e^{iq·ri}
    for (s, (ix, iy)) in enumerate(lat.sites)
        rx, ry = site_uc_pos(ix, iy)      # 用原胞 Bravais 位置，保证 n(q+G)=n(q)
        eiqr[s]   = cis(qx*rx + qy*ry)
        phase_sum += eiqr[s]
    end
    norm2 = sum(abs2(c) for c in values(coeffs))
    sq = 0.0
    for (F, c_F) in coeffs
        rho_q = zero(ComplexF64)
        tmp = F
        while tmp != 0
            s = trailing_zeros(tmp) + 1
            rho_q += eiqr[s]
            tmp &= tmp - 1
        end
        sq += abs2(c_F) * abs2(rho_q - n̄ * phase_sum)
    end
    return sq / (lat.Nuc * norm2)
end

# ── 均匀 BZ 网格扫描 n(q) ──
# N1, N2：沿 b1, b2 方向的采样数；n̄ = Np/Ns
# 返回 (qxs, qys, Sqs)
function structure_factor_bz(coeffs::Dict{Int64,ComplexF64},
                              lat::GenLat, N1::Int, N2::Int, n̄::Float64)
    b1 = collect(lat.b1); b2 = collect(lat.b2)
    Nq = N1 * N2
    qxs = zeros(Float64, Nq)
    qys = zeros(Float64, Nq)
    Sqs = zeros(Float64, Nq)
    idx = 0
    for h in 0:N1-1, k in 0:N2-1
        idx += 1
        qxs[idx] = h/N1 * b1[1] + k/N2 * b2[1]
        qys[idx] = h/N1 * b1[2] + k/N2 * b2[2]
        Sqs[idx] = structure_factor_at(coeffs, lat, qxs[idx], qys[idx], n̄)
    end
    return qxs, qys, Sqs
end

# ── 打印实空间密度报告 ──
function print_rspace_report(occ::Vector{Float64}, lat::GenLat, Np::Int)
    println("\n" * "="^60)
    println("实空间格点占据数 ⟨n_i⟩")
    n̄ = Np / lat.Ns
    @printf("  n̄ = %d/%d = %.6f\n", Np, lat.Ns, n̄)
    @printf("  ∑⟨n_i⟩ = %.6f (应为 %d)\n", sum(occ), Np)
    println("-"^60)
    println("  site   (ix, iy)    rx       ry      ⟨n_i⟩     δn=⟨n_i⟩-n̄")
    println("-"^60)
    for (s, (ix, iy)) in enumerate(lat.sites)
        rx, ry = site_phys_pos(ix, iy)
        @printf("  %3d   (%2d,%2d)  %6.3f  %6.3f   %.6f   %+.6f\n",
                s, ix, iy, rx, ry, occ[s], occ[s]-n̄)
    end
    occ_max, occ_min = maximum(occ), minimum(occ)
    @printf("\n  max ⟨n_i⟩ = %.6f (+%.4f σ)\n", occ_max, occ_max - n̄)
    @printf("  min ⟨n_i⟩ = %.6f (%+.4f σ)\n", occ_min, occ_min - n̄)
    @printf("  调制幅度 Δn = %.6f\n", occ_max - occ_min)
    println("="^60)
end
