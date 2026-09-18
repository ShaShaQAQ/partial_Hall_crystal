# ============================================================
# 实空间 Hopping 构造（通用版本，支持 GenLat）
# ============================================================
using LinearAlgebra

# k 空间哈密顿量（Eq.1 of PRL 113, 216404）
const A_PAPER = [[0.5, -sqrt(3)/2], [0.5, sqrt(3)/2], [-1.0, 0.0]]
const SIGMA_X = ComplexF64[0 1; 1 0]
const SIGMA_Y = ComplexF64[0 -im; im 0]
const SIGMA_Z = ComplexF64[1 0; 0 -1]
const IDENTITY_2 = Matrix{ComplexF64}(I, 2, 2)

function get_Hk_x_derivatives(k::Vector{Float64},
                               t1::Float64, t3::Float64)
    phases = [dot(k, a) for a in A_PAPER]
    g = [2t1 * cos(phase) for phase in phases]
    dg = [-2t1 * a[1] * sin(phase)
          for (a, phase) in zip(A_PAPER, phases)]
    d2g = [-2t1 * a[1]^2 * cos(phase)
           for (a, phase) in zip(A_PAPER, phases)]
    g0 = 2t3 * sum(cos(2phase) for phase in phases)
    dg0 = -4t3 * sum(a[1] * sin(2phase)
                     for (a, phase) in zip(A_PAPER, phases))
    d2g0 = -8t3 * sum(a[1]^2 * cos(2phase)
                      for (a, phase) in zip(A_PAPER, phases))
    assemble(v, v0) = (v[1] * SIGMA_X + v[2] * SIGMA_Y +
                       v[3] * SIGMA_Z + v0 * IDENTITY_2)
    return (Hk=assemble(g, g0),
            dHdkx=assemble(dg, dg0),
            d2Hdkx2=assemble(d2g, d2g0))
end

function get_Hk(k::Vector{Float64}, t1::Float64, t3::Float64)
    return get_Hk_x_derivatives(k, t1, t3).Hk
end

# 子格坐标由晶格约定给出。Tilted/legacy 几何使用负号 a2，
# RectLat4x6 使用与本地 ED 对齐的正号 a2。
sublat_pos(lat::GenLat) = [[0.0, 0.0], collect(lat.a2)]

# Fourier 变换得到实空间 hopping 矩阵 t(R)
# 使用 lat.kpoints（对不同超胞自动正确）
function fourier_to_real_from_bloch(lat::GenLat, bloch_matrix)
    Nk = length(lat.kpoints)
    a1_uc = collect(lat.a1)          # = (1,0)
    a2_uc = 2 .* collect(lat.a2)
    sublat = sublat_pos(lat)
    tR = Dict{NTuple{2,Int}, Matrix{ComplexF64}}()

    for (n1_uc, n2_uc) in lat.uc_trans
        R_cart = n1_uc .* a1_uc .+ n2_uc .* a2_uc
        tmat = zeros(ComplexF64, 2, 2)

        for k in lat.kpoints
            Hk = bloch_matrix(k)
            for α in 1:2, β in 1:2
                δ = sublat[α] .- sublat[β]
                phase = exp(-1im * dot(k, R_cart .+ δ))
                tmat[α, β] += phase * Hk[α, β]
            end
        end
        tmat ./= Nk
        tR[(n1_uc, n2_uc)] = tmat
    end
    return tR
end

function fourier_to_real(lat::GenLat, t1::Float64, t3::Float64)
    return fourier_to_real_from_bloch(lat, k -> get_Hk(k, t1, t3))
end

function fourier_to_real_Ax(lat::GenLat, t1::Float64,
                            t3::Float64, Ax::Float64)
    return fourier_to_real_from_bloch(
        lat, k -> get_Hk(k .+ [Ax, 0.0], t1, t3))
end

function fourier_to_real_x_derivatives(lat::GenLat,
                                       t1::Float64, t3::Float64)
    h = fourier_to_real_from_bloch(
        lat, k -> get_Hk_x_derivatives(k, t1, t3).Hk)
    jx = fourier_to_real_from_bloch(
        lat, k -> get_Hk_x_derivatives(k, t1, t3).dHdkx)
    kxx = fourier_to_real_from_bloch(
        lat, k -> get_Hk_x_derivatives(k, t1, t3).d2Hdkx2)
    return h, jx, kxx
end

function hopping_list_from_tR(lat::GenLat, tR;
                              threshold::Float64=1e-12)
    amplitudes = Dict{Tuple{Int,Int},ComplexF64}()
    for (source, (ix, iy)) in enumerate(lat.sites)
        source_orb = mod(iy, 2) + 1
        source_cell_y = div(iy - (source_orb - 1), 2)

        for ((n1_uc, n2_uc), tmat) in tR
            for target_orb in 1:2
                tix = ix + n1_uc
                tiy = 2 * (source_cell_y + n2_uc) + target_orb - 1
                target = lat.site_idx[canon_prim(lat, tix, tiy)]
                key = (target, source)
                amplitudes[key] = get(amplitudes, key, 0.0 + 0.0im) +
                                  tmat[target_orb, source_orb]
            end
        end
    end

    keys_sorted = sort!(collect(keys(amplitudes)))
    return [(target, source, amplitudes[(target, source)])
            for (target, source) in keys_sorted
            if abs(amplitudes[(target, source)]) >= threshold]
end

struct RealSpaceBond
    tgt::Int
    src::Int
    amp::ComplexF64
    winding_T2::Int
    displacement::NTuple{2,Float64}
end

function build_real_space_bonds(lat::GenLat, t1::Float64, t3::Float64)
    tR = fourier_to_real(lat, t1, t3)
    a1_uc = collect(lat.a1)
    a2_uc = 2 .* collect(lat.a2)
    sublat = sublat_pos(lat)
    bond_dict = Dict{Tuple{Int,Int},RealSpaceBond}()

    for (si, (ix, iy)) in enumerate(lat.sites)
        source_orb = mod(iy, 2) + 1
        source_cell_y = div(iy - (source_orb - 1), 2)

        for ((n1_uc, n2_uc), tmat) in tR
            Rcart = n1_uc .* a1_uc .+ n2_uc .* a2_uc
            for target_orb in 1:2
                amp = tmat[target_orb, source_orb]
                abs(amp) < 1e-12 && continue

                tix = ix + n1_uc
                tiy = 2 * (source_cell_y + n2_uc) + target_orb - 1

                target = canon_prim(lat, tix, tiy)
                haskey(lat.site_idx, target) || continue
                ti = lat.site_idx[target]

                haskey(bond_dict, (ti, si)) && continue

                delta = sublat[target_orb] .- sublat[source_orb]
                dcart = Rcart .+ delta
                d = (Float64(dcart[1]), Float64(dcart[2]))
                w2 = winding_T2(lat, tix, tiy)

                bond_dict[(ti, si)] = RealSpaceBond(ti, si, amp, w2, d)
                bond_dict[(si, ti)] = RealSpaceBond(
                    si, ti, conj(amp), -w2, (-d[1], -d[2]))
            end
        end
    end

    keys_sorted = sort!(collect(keys(bond_dict)))
    return [bond_dict[key] for key in keys_sorted]
end

# 构造完整实空间 hopping 列表（含 Hermitian 共轭）
# phi_y: 沿 T2 方向插入的磁通量（弧度）
function build_hops(lat::GenLat, t1::Float64, t3::Float64,
                    phi_y::Float64=0.0)
    return [(b.tgt, b.src, b.amp * cis(b.winding_T2 * phi_y))
            for b in build_real_space_bonds(lat, t1, t3)]
end


function build_hops_Ax(lat::GenLat, t1::Float64, t3::Float64,
                       Ax::Float64=0.0)
    tR = fourier_to_real_Ax(lat, t1, t3, Ax)
    return hopping_list_from_tR(lat, tR)
end


function build_hops_x_derivatives(lat::GenLat, t1::Float64, t3::Float64)
    h, jx, kxx = fourier_to_real_x_derivatives(lat, t1, t3)
    return (hopping_list_from_tR(lat, h),
            hopping_list_from_tR(lat, jx),
            hopping_list_from_tR(lat, kxx))
end
