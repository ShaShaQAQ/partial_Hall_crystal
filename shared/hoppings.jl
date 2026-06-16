# ============================================================
# 实空间 Hopping 构造（通用版本，支持 GenLat）
# ============================================================
using LinearAlgebra

# k 空间哈密顿量（Eq.1 of PRL 113, 216404）
function get_Hk(k::Vector{Float64}, t1::Float64, t3::Float64)
    σx = ComplexF64[0 1; 1 0]
    σy = ComplexF64[0 -im; im 0]
    σz = ComplexF64[1 0; 0 -1]
    I2 = Matrix{ComplexF64}(I, 2, 2)
    a_paper = [[0.5, -sqrt(3)/2], [0.5, sqrt(3)/2], [-1.0, 0.0]]
    g  = [2t1 * cos(dot(k, a)) for a in a_paper]
    g0 = 2t3 * sum(cos(2*dot(k, a)) for a in a_paper)
    return g[1]*σx + g[2]*σy + g[3]*σz + g0*I2
end

# 子格坐标：A=(0,0)，B=a2。这样 RectLat4x6 使用本地 ED 的
# a2=(1/2,sqrt(3)/2)，旧的 PHC 负号 a2 lattice 仍保持原坐标。
sublat_pos(lat::GenLat) = [[0.0, 0.0], collect(lat.a2)]

# Fourier 变换得到实空间 hopping 矩阵 t(R)
# 使用 lat.kpoints（对不同超胞自动正确）
function fourier_to_real(lat::GenLat, t1::Float64, t3::Float64)
    kpoints_list = lat.kpoints
    Nk = length(kpoints_list)

    a1_uc = collect(lat.a1)          # = (1,0)
    a2_uc = 2 .* collect(lat.a2)
    sublat = sublat_pos(lat)

    tR = Dict{NTuple{2,Int}, Matrix{ComplexF64}}()

    for (n1_uc, n2_uc) in lat.uc_trans
        R_cart = n1_uc .* a1_uc .+ n2_uc .* a2_uc
        tmat = zeros(ComplexF64, 2, 2)

        for k in kpoints_list
            Hk = get_Hk(k, t1, t3)
            for α in 1:2, β in 1:2
                δ = sublat[α] .- sublat[β]
                phase = exp(-1im * dot(k, R_cart .+ δ))
                tmat[α, β] += phase * Hk[α, β]
            end
        end
        tmat ./= Nk

        thresh = 1e-10
        for i in 1:2, j in 1:2
            abs(real(tmat[i,j])) < thresh && (tmat[i,j] -= real(tmat[i,j]))
            abs(imag(tmat[i,j])) < thresh && (tmat[i,j] -= 1im*imag(tmat[i,j]))
        end
        tR[(n1_uc, n2_uc)] = tmat
    end
    return tR
end

# 构造完整实空间 hopping 列表（含 Hermitian 共轭）
# phi_y: 沿 T2 方向插入的磁通量（弧度）
function build_hops(lat::GenLat, t1::Float64, t3::Float64,
                    phi_y::Float64=0.0)
    tR = fourier_to_real(lat, t1, t3)
    hop_dict = Dict{Tuple{Int,Int}, Tuple{ComplexF64,Int}}()

    for (si, (ix, iy)) in enumerate(lat.sites)
        sorb = (iy % 2 + 2) % 2

        for ((n1_uc, n2_uc), tmat) in tR
            for target_orb in 1:2
                amp = tmat[target_orb, sorb+1]
                abs(amp) < 1e-12 && continue

                target_orb_offset = target_orb - 1
                src_cell_y = div(iy - sorb, 2)
                dst_cell_x = ix + n1_uc
                dst_cell_y = src_cell_y + n2_uc
                tix = dst_cell_x
                tiy = 2*dst_cell_y + target_orb_offset

                r = canon_prim(lat, tix, tiy)   # 通用 canon
                haskey(lat.site_idx, r) || continue
                ti = lat.site_idx[r]

                haskey(hop_dict, (ti, si)) && continue

                w2 = winding_T2(lat, tix, tiy)  # 通用 winding
                hop_dict[(ti, si)] = (amp, w2)
                hop_dict[(si, ti)] = (conj(amp), -w2)
            end
        end
    end

    hops = Tuple{Int,Int,ComplexF64}[]
    for ((tgt, src), (amp, w2)) in hop_dict
        push!(hops, (tgt, src, amp * exp(1im * w2 * phi_y)))
    end
    return hops
end
