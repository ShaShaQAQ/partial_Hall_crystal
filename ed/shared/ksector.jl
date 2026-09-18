# ============================================================
# 动量扇区构造（通用版本）
#
# 相位直接查 lat.phase_table[m+1, ti]，无硬编码公式
# Bloch 态：|k,R> = (1/norm) Σ_T phase_table[m+1,T] * sgn_T * T|R>
# ============================================================

struct KSector
    m::Int
    reps::Vector{Int64}
    norms::Vector{Float64}
    orbit_data::Vector{Vector{Tuple{ComplexF64,Int64}}}
    fock2rep::Dict{Int64,Tuple{Int32,ComplexF64}}
end

function build_ksector(basis::Vector{Int64}, lat::GenLat, m::Int;
                       trans_maps::Union{Nothing,Matrix{Int}}=nothing)
    Ntrans = length(lat.Tnx)
    m_idx  = m + 1          # phase_table 的行索引（1-based）
    maps = trans_maps === nothing ? translation_site_maps(lat) : trans_maps
    size(maps) == (Ntrans, lat.Ns) ||
        throw(ArgumentError("translation map dimensions do not match lattice"))
    inversion_masks = translation_inversion_masks(maps, lat.Ns)
    visited = Set{Int64}()
    reps      = Int64[]
    norms     = Float64[]
    orb_datas = Vector{Vector{Tuple{ComplexF64,Int64}}}()
    f2r       = Dict{Int64,Tuple{Int32,ComplexF64}}()

    for F in basis
        F in visited && continue

        orb_coeffs = Dict{Int64,ComplexF64}()
        for ti in 1:Ntrans
            Fs, sgn = translate_fock_cached(
                F, maps, inversion_masks, ti, lat.Ns)
            ph = lat.phase_table[m_idx, ti] * sgn
            orb_coeffs[Fs] = get(orb_coeffs, Fs, 0.0+0im) + ph
        end

        rep = minimum(keys(orb_coeffs))
        if rep != F
            for Fs in keys(orb_coeffs); push!(visited, Fs); end
            continue
        end

        norm2 = sum(abs2, values(orb_coeffs))
        if norm2 < 1e-10
            for Fs in keys(orb_coeffs); push!(visited, Fs); end
            continue
        end

        norm = sqrt(norm2)
        rep_idx = Int32(length(reps) + 1)
        push!(reps, F)
        push!(norms, norm)

        orb_data_i = [(c_α, Fs) for (Fs, c_α) in orb_coeffs]
        push!(orb_datas, orb_data_i)

        for (Fs, c_α) in orb_coeffs
            f2r[Fs] = (rep_idx, conj(c_α) / norm)
        end
        for Fs in keys(orb_coeffs); push!(visited, Fs); end
    end

    return KSector(m, reps, norms, orb_datas, f2r)
end

function build_all_ksectors(basis::Vector{Int64}, lat::GenLat)
    secs = KSector[]
    maps = translation_site_maps(lat)
    for m in lat.ktab
        push!(secs, build_ksector(basis, lat, m; trans_maps=maps))
    end
    return secs
end
