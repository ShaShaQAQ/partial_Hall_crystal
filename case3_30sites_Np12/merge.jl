# ============================================================
# merge.jl — 汇总所有节点的 partial_*.jld2，计算结构因子，绘图
# 用法：julia merge.jl
# 前置条件：output/partial_0.jld2  partial_5.jld2  partial_10.jld2 均存在
# ============================================================
include("../shared/lattice.jl")
include("../shared/hoppings.jl")
include("../shared/basis.jl")
include("../shared/ksector.jl")
include("../shared/structure_factor.jl")

using Printf, Dates, JLD2

Np  = 12
V1  = 10.0
V2  = 2.0
V3  = 2.0
t1  = 1.0
t3  = 0.2

println("="^60)
println("Merge — 汇总 3 个节点的结果")
println("开始时间: ", now())
println("="^60)

t_total = time()

# ── 1. 重建晶格与完整扇区列表（供结构因子使用）──
print("[1/4] 重建晶格与扇区... "); flush(stdout)
t = @elapsed begin
    lat   = TiltedLat30()
    basis = gen_basis(lat.Ns, Np)
    secs  = KSector[]
    for m in lat.ktab
        push!(secs, build_ksector(basis, lat, m))
    end
end
@printf("done  Ns=%d  Nuc=%d  (%.2f s)\n", lat.Ns, lat.Nuc, t)

# ── 2. 加载所有 partial 文件 ──
println("[2/4] 加载 partial JLD2 文件...")
partial_files = [
    ("output/partial_0.jld2",  "k=0..4"),
    ("output/partial_5.jld2",  "k=5..9"),
    ("output/partial_10.jld2", "k=10..14"),
]

all_ev   = Tuple{Int,Float64}[]
gs_vecs  = Dict{Int,Vector{ComplexF64}}()

for (fpath, label) in partial_files
    isfile(fpath) || error("找不到文件: $fpath")
    data = load(fpath)
    append!(all_ev, data["ev_pairs"])
    merge!(gs_vecs, data["gs_vecs"])
    @printf("  ✓ %s  (%s)\n", fpath, label)
end
sort!(all_ev, by=x->x[2])
E0   = all_ev[1][2]
gs_k = all_ev[1][1]
@printf("  全局基态: k=%d  E0=%.10f\n", gs_k, E0)

# ── 3. 输出能谱文件 ──
println("[3/4] 输出结果文件...")
mkpath("output")

open("output/spectrum_Np$(Np).dat","w") do f
    println(f,"# k  E-E0  [30sites Np=$Np V1=$V1 V2=$V2 V3=$V3 t'=$t3 multi-node]")
    for (m,e) in all_ev
        @printf(f,"%d  %.10f\n", m, e-E0)
    end
end
println("  能谱保存: output/spectrum_Np$(Np).dat")

# ── 结构因子 ──
gs_sec = secs[findfirst(s -> s.m == gs_k, secs)]
print("  全局基态结构因子... "); flush(stdout)
t_sq = @elapsed begin
    gs_vec = gs_vecs[gs_k]
    results_ord, results_srt, norm2, Np_check =
        structure_factor_survey(gs_sec, gs_vec, lat, Np)
end
@printf("done (%.1f s)\n", t_sq)
print_sq_report(results_srt, norm2, Np_check, Np, lat.Ns, gs_k)

open("output/sq_Np$(Np).dat","w") do f
    println(f,"# ik  N(k)  kx  ky  [30sites Np=$Np V1=$V1 t'=$t3]")
    println(f,"# 基态扇区 k=$gs_k  E0=$(round(E0,digits=8))")
    for (m, sq, kx, ky) in results_ord
        @printf(f,"%d  %.6f  %.6f  %.6f\n", m, sq, kx, ky)
    end
end
println("  结构因子保存: output/sq_Np$(Np).dat")

print("  全扇区结构因子矩阵... "); flush(stdout)
t_all = @elapsed open("output/sq_all_Np$(Np).dat","w") do f
    println(f,"# sector_k  q  N(q)  kx  ky")
    for sec in secs
        haskey(gs_vecs, sec.m) || continue
        res_ord, _, _, _ = structure_factor_survey(sec, gs_vecs[sec.m], lat, Np)
        for (m, sq, kx, ky) in res_ord
            @printf(f,"%d  %d  %.6f  %.6f  %.6f\n", sec.m, m, sq, kx, ky)
        end
    end
end
@printf("done (%.1f s)\n", t_all)
println("  全扇区结构因子保存: output/sq_all_Np$(Np).dat")

# ── 4. 绘图 ──
println("[4/4] 绘图...")
# plot.jl 读取当前目录下的 dat 文件，切换到 output/ 再 include
cd("output")
include("../plot.jl")
cd("..")

@printf("\n总耗时: %.1f s (%.2f min)\n", time()-t_total, (time()-t_total)/60)
println("完成时间: ", now())
println("完成!")
