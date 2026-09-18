# ============================================================
# validate_tilted30.jl
# TiltedLat30 投影 ED vs 全 2-band ED 对比
#   Np=12  ν = 12/15 = 4/5（空穴填充 1/5）
#   Np=13  ν = 13/15
#   参数：t=1, t'=0.2, V1=10, V2=2, V3=2
# 全 ED 参考：case3_30sites_Np12/output_Np{12,13}/
# ============================================================
include("../../shared/lattice.jl")
include("../../shared/hoppings.jl")
include("../shared_proj/band.jl")
include("../shared_proj/proj_basis.jl")
include("../shared_proj/proj_interaction.jl")
include("../shared_proj/proj_hamiltonian.jl")
include("../shared_proj/proj_structure_factor.jl")

using Printf, LinearAlgebra

# ── 全局参数 ──
t1  = 1.0
t3  = 0.2
V1  = 10.0
V2  = 2.0
V3  = 2.0

# ── 构建晶格（仅一次）──
println("="^60)
print("构建 TiltedLat30... "); flush(stdout)
lat  = TiltedLat30()
Nk   = lat.Nuc
eps, uvec = compute_band(lat, t1, t3)    # 最低 Chern 带本征矢与单粒子能量

@printf("done  Nk=%d  Ns=%d\n", Nk, lat.Ns)
@printf("带宽: [%.6f, %.6f]\n", minimum(eps), maximum(eps))

# 预计算 k_add 与 W（共用于所有 Np）
print("预计算 k_add... "); flush(stdout)
kadd = build_k_add(lat)
println("done")

print("预计算 W (V1=$V1, V2=$V2, V3=$V3)... "); flush(stdout)
W = build_Vq_V123(lat, uvec, V1, V2, V3)
println("done")

println("="^60)

# ============================================================
function run_np(Np::Int, full_ed_dir::String)
    label = "Np=$Np  ν=$(Np)/$Nk"
    println("\n" * "="^60)
    println("投影 ED: $label")
    println("="^60)

    # ── 基矢与扇区 ──
    sectors = build_proj_sectors(Nk, Np)
    total   = sum(length(v) for v in values(sectors))
    @printf("Hilbert 空间: C(%d,%d) = %d  (15 个扇区)\n", Nk, Np, total)
    for M in 0:Nk-1
        @printf("  M=%2d: dim=%d\n", M, length(sectors[M]))
    end

    # ── 投影 ED 求解 ──
    evp, all_vecs = solve_proj_ed(sectors, eps, W, kadd, Nk)
    E0p = evp[1][2]
    @printf("\n投影 ED: E0=%.6f  基态扇区 M=%d\n", E0p, evp[1][1])

    # 打印前20个能级（相对基态）
    println("\n前 20 个能级（ΔE 相对基态）:")
    @printf("  %-4s  %-6s  %-12s\n", "rank", "k", "ΔE_proj")
    println("  " * "-"^24)
    for (i, (m, e)) in enumerate(evp[1:min(20,end)])
        @printf("  %-4d  %-6d  %.10f\n", i, m, e - E0p)
    end

    # ── 全 ED 能谱读取与对比 ──
    full_path = joinpath(full_ed_dir, "spectrum_Np$(Np).dat")
    if isfile(full_path)
        evf = Tuple{Int,Float64}[]
        for l in readlines(full_path)
            startswith(l,"#") && continue
            p = split(strip(l)); length(p)<2 && continue
            push!(evf, (parse(Int,p[1]), parse(Float64,p[2])))
        end
        E0f = 0.0  # 已经是 E-E0 格式

        println("\n投影 ED vs 全 2-band ED（前20态）:")
        @printf("  %-4s  %-6s  %-14s  %-14s  %-12s\n",
            "rank", "k_proj", "ΔE_proj", "ΔE_full", "偏差")
        println("  " * "-"^56)
        n_show = min(20, length(evp), length(evf))
        for i in 1:n_show
            m_p, e_p = evp[i]
            m_f, e_f = evf[i]
            @printf("  %-4d  %-6d  %.10f  %.10f  %+.4e\n",
                i, m_p, e_p - E0p, e_f, (e_p-E0p) - e_f)
        end

        # GS 简并检验
        gs_k_p = [m for (m,e) in evp if abs(e - E0p) < 1e-8]
        gs_k_f = [m for (m,e) in evf if e < 1e-8]
        @printf("\n投影 ED 基态扇区: %s\n", join(gs_k_p, ", "))
        @printf("全 2-band 基态扇区: %s\n", join(gs_k_f, ", "))

        # spread 与 gap
        ep_sorted = sort(evp[1:end], by=x->x[2])
        n_gs = length(gs_k_p)
        if length(ep_sorted) > n_gs
            spread_p = ep_sorted[n_gs][2]   - ep_sorted[1][2]
            gap_p    = ep_sorted[n_gs+1][2] - ep_sorted[n_gs][2]
            @printf("投影 ED: GS散布=%.6f  gap=%.6f\n", spread_p, gap_p)
        end
    else
        println("未找到全 ED 参考: $full_path")
        evf = Tuple{Int,Float64}[]
    end

    # ── 保存投影能谱 ──
    out_spec = "spectrum_proj_tilted30_Np$(Np).dat"
    open(out_spec, "w") do f
        println(f, "# k  E-E0  [proj ED TiltedLat30 Np=$Np V1=$V1 V2=$V2 V3=$V3 t3=$t3]")
        for (m, e) in evp
            @printf(f, "%d  %.10f\n", m, e - E0p)
        end
    end
    println("\n保存能谱: $out_spec")

    # ── 结构因子 ──
    M_gs   = evp[1][1]
    gs_vec = complex(all_vecs[M_gs][:, 1])
    sq_p   = compute_proj_sq(gs_vec, sectors[M_gs], kadd, Nk, uvec)

    out_sq = "sq_proj_tilted30_Np$(Np).dat"
    open(out_sq, "w") do f
        println(f, "# ik  N_proj(k)  kx  ky  [proj ED TiltedLat30 Np=$Np V1=$V1 V2=$V2 V3=$V3 GS_sector=$M_gs]")
        for qi in 0:Nk-1
            kp = lat.kpoints[qi+1]
            @printf(f, "%d  %.8f  %.6f  %.6f\n", qi, sq_p[qi+1], kp[1], kp[2])
        end
    end

    # 打印结构因子前5强峰
    sq_sorted = sort(collect(enumerate(sq_p)), by=x->-x[2])
    println("N_proj(q) 前5强峰：")
    for (qi1, val) in sq_sorted[1:min(5,end)]
        kp = lat.kpoints[qi1]
        @printf("  q=%2d  N=%.5f  (kx=%.4f, ky=%.4f)\n", qi1-1, val, kp[1], kp[2])
    end
    println("保存结构因子: $out_sq")

    # 全扇区结构因子矩阵
    sq_all = compute_proj_sq_allsectors(all_vecs, sectors, kadd, Nk, uvec)
    out_sq_all = "sq_all_proj_tilted30_Np$(Np).dat"
    open(out_sq_all, "w") do f
        println(f, "# sector_k  q  N_proj(q)  [proj ED TiltedLat30 Np=$Np V1=$V1 V2=$V2 V3=$V3]")
        for M in 0:Nk-1
            for qi in 0:Nk-1
                @printf(f, "%d  %d  %.8f\n", M, qi, sq_all[M+1, qi+1])
            end
        end
    end
    println("保存全扇区结构因子: $out_sq_all")

    # 读取全 ED 结构因子（如有）
    sq_full_path = joinpath(full_ed_dir, "sq_Np$(Np).dat")
    if isfile(sq_full_path)
        sq_full = Dict{Int,Float64}()
        for l in readlines(sq_full_path)
            startswith(l,"#") && continue
            p = split(strip(l)); length(p)<2 && continue
            sq_full[parse(Int,p[1])] = parse(Float64,p[2])
        end
        sf_sorted = sort(collect(sq_full), by=x->-x[2])
        println("全 2-band ED N(q) 前5强峰（旧 kpoints 仅供参考）：")
        for (q, v) in sf_sorted[1:min(5,end)]
            @printf("  q=%2d  N_full=%.5f\n", q, v)
        end
    end

    return evp, all_vecs, sq_p
end
# ============================================================

full_dir12 = "../../cases/case3_30sites_Np12/output_Np12"
full_dir13 = "../../cases/case3_30sites_Np12/output_Np13"

evp12, vecs12, sq12 = run_np(12, full_dir12)
evp13, vecs13, sq13 = run_np(13, full_dir13)

println("\n" * "="^60)
println("全部完成。")
println("输出文件：")
println("  spectrum_proj_tilted30_Np{12,13}.dat")
println("  sq_proj_tilted30_Np{12,13}.dat")
println("  sq_all_proj_tilted30_Np{12,13}.dat")
