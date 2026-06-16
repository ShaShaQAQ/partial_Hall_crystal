# ============================================================
# 统一晶格模块 GenLat
# 支持：
#   TiltedLat30() — 倾斜 (4×4-1) 团簇，30 sites，15 unit cells
#   RectLat3x6()  — 3×6 矩形超胞，18 sites，9 unit cells
#
# 关键：phase_table[m+1, ti] 预计算好，ksector 直接查表
# ============================================================

struct GenLat
    Ns::Int                                    # 格点总数
    Nuc::Int                                   # 单元格数 = k 点数
    a1::NTuple{2,Float64}                      # 三角晶格格矢（物理坐标）
    a2::NTuple{2,Float64}
    b1::NTuple{2,Float64}                      # 倒格矢
    b2::NTuple{2,Float64}
    sites::Vector{NTuple{2,Int}}               # site_index → (ix,iy)
    site_idx::Dict{NTuple{2,Int},Int}          # (ix,iy) → site_index
    uc_trans::Vector{NTuple{2,Int}}            # 单元格平移 (n1_uc, n2_uc)
    Tnx::Vector{Int}                           # 三角坐标 x 平移
    Tny::Vector{Int}                           # 三角坐标 y 平移（= 2*n2_uc）
    ktab::Vector{Int}                          # k 标签 0:Nuc-1
    phase_table::Matrix{ComplexF64}            # [Nuc × Ntrans]: phase_table[m+1, ti]
    kpoints::Vector{Vector{Float64}}           # Nuc 个物理 k 矢（用于 Fourier 变换）
    T1::NTuple{2,Int}                          # 超胞矢量 T1（三角坐标）
    T2::NTuple{2,Int}                          # 超胞矢量 T2（三角坐标）
    nb1::Matrix{Int}                           # 最近邻 Ns×6
    nb2::Matrix{Int}                           # 次近邻 Ns×6
    nb3::Matrix{Int}                           # 第三近邻 Ns×6
end

# ── 通用 canon_prim: 将 (ix,iy) 规约到超胞基本域 ──
# 超胞矩阵 M = [T1|T2]，det(M) = Ns
# 逆 * det: [[T2y,-T2x],[-T1y,T1x]]
@inline function canon_prim(lat::GenLat, ix::Int, iy::Int)
    T1x, T1y = lat.T1
    T2x, T2y = lat.T2
    n1 = fld(T2y*ix - T2x*iy, lat.Ns)
    n2 = fld(-T1y*ix + T1x*iy, lat.Ns)
    return (ix - T1x*n1 - T2x*n2, iy - T1y*n1 - T2y*n2)
end

# ── 沿 T2 方向越界次数（谱流通量用）──
@inline function winding_T2(lat::GenLat, tix::Int, tiy::Int)
    T1x, T1y = lat.T1
    return fld(-T1y*tix + T1x*tiy, lat.Ns)
end

# ============================================================
# 构造函数 1：倾斜 (4×4-1) 30 格点团簇
# T1_tri=(4,2)，T2_tri=(1,8)，det=30，Nuc=15
# 相位公式：exp(-2πi/15 * m * (4*n1_uc - n2_uc))
# ============================================================
function TiltedLat30()
    T1 = (4, 2);  T2 = (1, 8);  Ns = 30;  Nuc = 15

    # 局部 canon（构造时用）
    function cp(ix, iy)
        n1 = fld(8ix - iy, 30)
        n2 = fld(-2ix + 4iy, 30)
        return (ix - 4n1 - n2, iy - 2n1 - 8n2)
    end

    # 枚举 30 个格点
    seen = Dict{NTuple{2,Int},Bool}()
    for ix in -6:16, iy in -4:22
        seen[cp(ix,iy)] = true
    end
    sites = sort(collect(keys(seen)))
    @assert length(sites) == Ns "格点数错误: $(length(sites))"
    site_idx = Dict(s => i for (i,s) in enumerate(sites))

    # 物理格矢与倒格矢
    a1 = (1.0, 0.0);  a2 = (-0.5, sqrt(3)/2)
    area = a1[1]*a2[2] - a1[2]*a2[1]
    b1 = (2π*a2[2]/area, -2π*a2[1]/area)
    b2 = (-2π*a1[2]/area,  2π*a1[1]/area)

    # 15 个单元格平移（canon_uc：T1_uc=(4,1), T2_uc=(1,4)）
    function cu(p, q)
        n1 = fld(4p - q, 15);  n2 = fld(-p + 4q, 15)
        return (p - 4n1 - n2, q - n1 - 4n2)
    end
    uc_seen = Dict{NTuple{2,Int},Bool}()
    for p in -8:12, q in -8:12;  uc_seen[cu(p,q)] = true;  end
    uc_trans = sort(collect(keys(uc_seen)))
    @assert length(uc_trans) == 15
    Tnx = [t[1] for t in uc_trans]
    Tny = [2*t[2] for t in uc_trans]

    # 相位表
    phase_table = Matrix{ComplexF64}(undef, Nuc, Nuc)
    for (mi, m) in enumerate(0:14)
        for ti in 1:Nuc
            n1_uc = Tnx[ti];  n2_uc = Tny[ti] ÷ 2
            phase_table[mi, ti] = cis(-2π/15 * m * (4*n1_uc - n2_uc))
        end
    end

    # k 点（用于 Fourier 变换）
    b1c = collect(b1);  b2c = collect(b2)
    G1_sc = (4 .* b1c .- b2c ./ 2) ./ 15
    kpoints = [m .* G1_sc for m in 0:14]

    # 近邻矩阵
    d1 = [(1,0),(-1,0),(1,1),(-1,-1),(0,1),(0,-1)]
    d2 = [(2,1),(-2,-1),(1,2),(-1,-2),(1,-1),(-1,1)]
    d3 = [(2,0),(-2,0),(0,2),(0,-2),(2,2),(-2,-2)]
    function make_nb(ds)
        M = Matrix{Int}(undef, Ns, 6)
        for (si,(ix,iy)) in enumerate(sites)
            for (k,(dx,dy)) in enumerate(ds)
                M[si,k] = site_idx[cp(ix+dx, iy+dy)]
            end
        end
        M
    end

    return GenLat(Ns, Nuc, a1, a2, b1, b2, sites, site_idx,
                  uc_trans, Tnx, Tny, collect(0:14), phase_table, kpoints,
                  T1, T2, make_nb(d1), make_nb(d2), make_nb(d3))
end

# ============================================================
# 构造函数 2：3×6 矩形超胞，18 格点，9 单元格
# T1_tri=(3,0)，T2_tri=(0,6)，det=18，Nuc=9
# k 标签 m = m1 + 3*m2，m1,m2 ∈ {0,1,2}
# 相位公式：exp(-2πi/3 * ((m%3)*n1_uc + (m÷3)*n2_uc))
# ============================================================
function RectLat3x6()
    T1 = (3, 0);  T2 = (0, 6);  Ns = 18;  Nuc = 9

    # 局部 canon
    function cp(ix, iy)
        n1 = fld(ix, 3);  n2 = fld(iy, 6)
        return (ix - 3n1, iy - 6n2)
    end

    # 枚举 18 个格点
    seen = Dict{NTuple{2,Int},Bool}()
    for ix in -4:8, iy in -6:12
        seen[cp(ix,iy)] = true
    end
    sites = sort(collect(keys(seen)))
    @assert length(sites) == Ns "格点数错误: $(length(sites))"
    site_idx = Dict(s => i for (i,s) in enumerate(sites))

    # 物理格矢与倒格矢（与 TiltedLat30 相同）
    a1 = (1.0, 0.0);  a2 = (-0.5, sqrt(3)/2)
    area = a1[1]*a2[2] - a1[2]*a2[1]
    b1 = (2π*a2[2]/area, -2π*a2[1]/area)
    b2 = (-2π*a1[2]/area,  2π*a1[1]/area)

    # 9 个单元格平移（T1_uc=(3,0), T2_uc=(0,3)）
    function cu(p, q)
        n1 = fld(p, 3);  n2 = fld(q, 3)
        return (p - 3n1, q - 3n2)
    end
    uc_seen = Dict{NTuple{2,Int},Bool}()
    for p in -4:6, q in -4:6;  uc_seen[cu(p,q)] = true;  end
    uc_trans = sort(collect(keys(uc_seen)))
    @assert length(uc_trans) == 9 "单元格平移数错误: $(length(uc_trans))"
    Tnx = [t[1] for t in uc_trans]
    Tny = [2*t[2] for t in uc_trans]

    # 相位表：m = m1 + 3*m2
    phase_table = Matrix{ComplexF64}(undef, Nuc, Nuc)
    for (mi, m) in enumerate(0:8)
        m1 = m % 3;  m2 = m ÷ 3
        for ti in 1:Nuc
            n1_uc = Tnx[ti];  n2_uc = Tny[ti] ÷ 2
            phase_table[mi, ti] = cis(-2π/3 * (m1*n1_uc + m2*n2_uc))
        end
    end

    # k 点（Fourier 变换用）
    # G1_sc = B1_uc/3 = b1/3，G2_sc = B2_uc/3 = b2/6
    b1c = collect(b1);  b2c = collect(b2)
    G1_sc = b1c ./ 3
    G2_sc = b2c ./ 6
    kpoints = [(m%3) .* G1_sc .+ (m÷3) .* G2_sc for m in 0:8]

    # 近邻矩阵（同位移向量）
    d1 = [(1,0),(-1,0),(1,1),(-1,-1),(0,1),(0,-1)]
    d2 = [(2,1),(-2,-1),(1,2),(-1,-2),(1,-1),(-1,1)]
    d3 = [(2,0),(-2,0),(0,2),(0,-2),(2,2),(-2,-2)]
    function make_nb(ds)
        M = Matrix{Int}(undef, Ns, 6)
        for (si,(ix,iy)) in enumerate(sites)
            for (k,(dx,dy)) in enumerate(ds)
                M[si,k] = site_idx[cp(ix+dx, iy+dy)]
            end
        end
        M
    end

    return GenLat(Ns, Nuc, a1, a2, b1, b2, sites, site_idx,
                  uc_trans, Tnx, Tny, collect(0:8), phase_table, kpoints,
                  T1, T2, make_nb(d1), make_nb(d2), make_nb(d3))
end

# ============================================================
# 构造函数 3：5×6 sites（5×3 原胞），30 格点，15 单元格
# T1_tri=(5,0)，T2_tri=(0,6)，det=30，Nuc=15
# k 标签 m = m1 + 5*m2，m1∈{0..4}，m2∈{0..2}
# 相位公式：exp(-2πi * (m1*n1_uc/5 + m2*n2_uc/3))
# 填充 Np=12 → ν = 12/15 = 4/5
# ============================================================
function RectLat5x6()
    T1 = (5, 0);  T2 = (0, 6);  Ns = 30;  Nuc = 15

    function cp(ix, iy)
        n1 = fld(ix, 5);  n2 = fld(iy, 6)
        return (ix - 5n1, iy - 6n2)
    end

    seen = Dict{NTuple{2,Int},Bool}()
    for ix in -6:12, iy in -6:12
        seen[cp(ix,iy)] = true
    end
    sites = sort(collect(keys(seen)))
    @assert length(sites) == Ns "格点数错误: $(length(sites))"
    site_idx = Dict(s => i for (i,s) in enumerate(sites))

    a1 = (1.0, 0.0);  a2 = (-0.5, sqrt(3)/2)
    area = a1[1]*a2[2] - a1[2]*a2[1]
    b1 = (2π*a2[2]/area, -2π*a2[1]/area)
    b2 = (-2π*a1[2]/area,  2π*a1[1]/area)

    # T1_uc=(5,0)，T2_uc=(0,3) → 15 个单元格平移
    function cu(p, q)
        n1 = fld(p, 5);  n2 = fld(q, 3)
        return (p - 5n1, q - 3n2)
    end
    uc_seen = Dict{NTuple{2,Int},Bool}()
    for p in -6:12, q in -4:6;  uc_seen[cu(p,q)] = true;  end
    uc_trans = sort(collect(keys(uc_seen)))
    @assert length(uc_trans) == 15 "单元格平移数错误: $(length(uc_trans))"
    Tnx = [t[1] for t in uc_trans]
    Tny = [2*t[2] for t in uc_trans]

    # 相位表：m = m1 + 5*m2
    phase_table = Matrix{ComplexF64}(undef, Nuc, Nuc)
    for (mi, m) in enumerate(0:14)
        m1 = m % 5;  m2 = m ÷ 5
        for ti in 1:Nuc
            n1_uc = Tnx[ti];  n2_uc = Tny[ti] ÷ 2
            phase_table[mi, ti] = cis(-2π * (m1*n1_uc/5 + m2*n2_uc/3))
        end
    end

    # k 点：G1_sc = b1/5，G2_sc = b2/6
    b1c = collect(b1);  b2c = collect(b2)
    G1_sc = b1c ./ 5
    G2_sc = b2c ./ 6
    kpoints = [(m%5) .* G1_sc .+ (m÷5) .* G2_sc for m in 0:14]

    d1 = [(1,0),(-1,0),(1,1),(-1,-1),(0,1),(0,-1)]
    d2 = [(2,1),(-2,-1),(1,2),(-1,-2),(1,-1),(-1,1)]
    d3 = [(2,0),(-2,0),(0,2),(0,-2),(2,2),(-2,-2)]
    function make_nb(ds)
        M = Matrix{Int}(undef, Ns, 6)
        for (si,(ix,iy)) in enumerate(sites)
            for (k,(dx,dy)) in enumerate(ds)
                M[si,k] = site_idx[cp(ix+dx, iy+dy)]
            end
        end
        M
    end

    return GenLat(Ns, Nuc, a1, a2, b1, b2, sites, site_idx,
                  uc_trans, Tnx, Tny, collect(0:14), phase_table, kpoints,
                  T1, T2, make_nb(d1), make_nb(d2), make_nb(d3))
end

# ============================================================
# 构造函数 4：4×6 sites（4×3 原胞），24 格点，12 单元格
# T1_tri=(4,0)，T2_tri=(0,6)，det=24，Nuc=12
# k 标签 m = m1 + 4*m2，m1∈{0..3}，m2∈{0..2}
# 相位公式：exp(-2πi * (m1*n1_uc/4 + m2*n2_uc/3))
# 填充 Np=4 → ν = 4/12 = 1/3（下 Chern 带）
# ============================================================
function RectLat4x6()
    T1 = (4, 0);  T2 = (0, 6);  Ns = 24;  Nuc = 12

    function cp(ix, iy)
        n1 = fld(ix, 4);  n2 = fld(iy, 6)
        return (ix - 4n1, iy - 6n2)
    end

    seen = Dict{NTuple{2,Int},Bool}()
    for ix in -5:10, iy in -6:12
        seen[cp(ix,iy)] = true
    end
    sites = sort(collect(keys(seen)))
    @assert length(sites) == Ns "格点数错误: $(length(sites))"
    site_idx = Dict(s => i for (i,s) in enumerate(sites))

    # 与本地 ED 保持同一套三角晶格基矢约定：
    # a1=(1,0), a2=(1/2,sqrt(3)/2)。
    a1 = (1.0, 0.0);  a2 = (0.5, sqrt(3)/2)
    area = a1[1]*a2[2] - a1[2]*a2[1]
    b1 = (2π*a2[2]/area, -2π*a2[1]/area)
    b2 = (-2π*a1[2]/area,  2π*a1[1]/area)

    # T1_uc=(4,0)，T2_uc=(0,3) → 12 个单元格平移
    function cu(p, q)
        n1 = fld(p, 4);  n2 = fld(q, 3)
        return (p - 4n1, q - 3n2)
    end
    uc_seen = Dict{NTuple{2,Int},Bool}()
    for p in -5:10, q in -4:6;  uc_seen[cu(p,q)] = true;  end
    uc_trans = sort(collect(keys(uc_seen)))
    @assert length(uc_trans) == 12 "单元格平移数错误: $(length(uc_trans))"
    Tnx = [t[1] for t in uc_trans]
    Tny = [2*t[2] for t in uc_trans]

    # 相位表：m = m1 + 4*m2
    phase_table = Matrix{ComplexF64}(undef, Nuc, Nuc)
    for (mi, m) in enumerate(0:11)
        m1 = m % 4;  m2 = m ÷ 4
        for ti in 1:Nuc
            n1_uc = Tnx[ti];  n2_uc = Tny[ti] ÷ 2
            phase_table[mi, ti] = cis(-2π * (m1*n1_uc/4 + m2*n2_uc/3))
        end
    end

    # k 点：G1_sc = b1/4，G2_sc = b2/6
    b1c = collect(b1);  b2c = collect(b2)
    G1_sc = b1c ./ 4
    G2_sc = b2c ./ 6
    kpoints = [(m%4) .* G1_sc .+ (m÷4) .* G2_sc for m in 0:11]

    d1 = [(1,0),(-1,0),(0,1),(0,-1),(-1,1),(1,-1)]
    d2 = [(1,1),(-1,-1),(-1,2),(1,-2),(-2,1),(2,-1)]
    d3 = [(2,0),(-2,0),(0,2),(0,-2),(-2,2),(2,-2)]
    function make_nb(ds)
        M = Matrix{Int}(undef, Ns, 6)
        for (si,(ix,iy)) in enumerate(sites)
            for (k,(dx,dy)) in enumerate(ds)
                M[si,k] = site_idx[cp(ix+dx, iy+dy)]
            end
        end
        M
    end

    return GenLat(Ns, Nuc, a1, a2, b1, b2, sites, site_idx,
                  uc_trans, Tnx, Tny, collect(0:11), phase_table, kpoints,
                  T1, T2, make_nb(d1), make_nb(d2), make_nb(d3))
end
