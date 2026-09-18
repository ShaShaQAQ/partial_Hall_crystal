# ============================================================
# Fock 态基矢生成与平移操作（含费米符号）
# ============================================================

# 生成所有 C(Ns, Np) 个 Fock 态（按升序排列）
function gen_basis(Ns::Int, Np::Int)
    basis = Int64[]
    function rec!(s, start, left)
        left == 0 && (push!(basis, s); return)
        for i in start:Ns-left+1
            rec!(s | (Int64(1)<<(i-1)), i+1, left-1)
        end
    end
    rec!(Int64(0), 1, Np)
    sort!(basis)
    return basis
end

# 应用 c†_i c_j：从 site j 移到 site i
@inline function hop(state::Int64, i::Int, j::Int)
    (state>>(j-1))&1 == 0 && return (Int64(0), 0)
    (state>>(i-1))&1 == 1 && return (Int64(0), 0)
    new_state = state ⊻ (Int64(1)<<(j-1)) ⊻ (Int64(1)<<(i-1))
    lo, hi = min(i,j), max(i,j)
    between_mask = 0
    if hi > lo + 1
        between_mask = ((Int64(1)<<(hi-1)) - 1) ⊻ ((Int64(1)<<lo) - 1)
    end
    sgn = isodd(count_ones(state & between_mask)) ? -1 : 1
    return (new_state, sgn)
end

# 将 Fock 态沿 (nx, ny) 平移，返回新态和费米符号
function translate_fock(state::Int64, lat::GenLat, nx::Int, ny::Int)
    occ = Int[]
    for s in 1:lat.Ns
        (state>>(s-1))&1 == 1 && push!(occ, s)
    end

    new_occ = Int[]
    for s in occ
        ix, iy = lat.sites[s]
        r = canon_prim(lat, ix+nx, iy+ny)   # 使用通用 canon
        push!(new_occ, lat.site_idx[r])
    end

    perm = sortperm(new_occ)
    n_inv = 0
    np = length(perm)
    for i in 1:np-1, j in i+1:np
        n_inv += (perm[i] > perm[j])
    end
    sgn = isodd(n_inv) ? -1 : 1

    new_state = Int64(0)
    for s in new_occ
        new_state |= Int64(1)<<(s-1)
    end
    return new_state, sgn
end

# Precompute the site permutation for every unit-cell translation.
function translation_site_maps(lat::GenLat)
    maps = Matrix{Int}(undef, length(lat.uc_trans), lat.Ns)
    for (ti, (nx, ny)) in enumerate(zip(lat.Tnx, lat.Tny))
        @inbounds for s in 1:lat.Ns
            ix, iy = lat.sites[s]
            maps[ti, s] = lat.site_idx[canon_prim(lat, ix + nx, iy + ny)]
        end
    end
    return maps
end

# Cache the translated-order inversions which determine the fermion sign.
function translation_inversion_masks(maps::Matrix{Int}, Ns::Int)
    masks = Matrix{UInt64}(undef, size(maps, 1), Ns)
    @inbounds for ti in axes(maps, 1), s in 1:Ns
        translated_s = maps[ti, s]
        mask = UInt64(0)
        for t in (s + 1):Ns
            maps[ti, t] < translated_s &&
                (mask |= UInt64(1) << (t - 1))
        end
        masks[ti, s] = mask
    end
    return masks
end

"""Allocation-free Fock translation using precomputed site maps and signs."""
@inline function translate_fock_cached(
        state::Int64,
        maps::Matrix{Int},
        inversion_masks::Matrix{UInt64},
        translation_index::Int,
        Ns::Int)
    new_state = Int64(0)
    occupied = state
    while occupied != 0
        s = trailing_zeros(occupied) + 1
        new_state |= Int64(1) << (maps[translation_index, s] - 1)
        occupied &= occupied - 1
    end

    parity = false
    occupied = state
    state_unsigned = UInt64(state)
    while occupied != 0
        s = trailing_zeros(occupied) + 1
        parity = xor(
            parity,
            isodd(count_ones(
                state_unsigned & inversion_masks[translation_index, s])))
        occupied &= occupied - 1
    end
    return new_state, parity ? -1 : 1
end
