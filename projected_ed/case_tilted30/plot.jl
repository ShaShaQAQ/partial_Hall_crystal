# ============================================================
# plot.jl — TiltedLat30 投影 ED vs 全 ED 能谱 & 结构因子图
# 用法：julia plot.jl（需先跑完 validate_tilted30.jl）
# ============================================================
using Plots, Printf
gr()

Nk  = 15
klabels = string.(0:Nk-1)

# ── 工具函数 ──
function read_spectrum(path)
    ks=Int[]; es=Float64[]
    for l in readlines(path)
        startswith(l,"#") && continue
        p=split(strip(l)); length(p)<2 && continue
        push!(ks,parse(Int,p[1])); push!(es,parse(Float64,p[2]))
    end
    return ks, es
end

function read_sq(path)
    qs=Int[]; ns=Float64[]
    for l in readlines(path)
        startswith(l,"#") && continue
        p=split(strip(l)); length(p)<2 && continue
        push!(qs,parse(Int,p[1])); push!(ns,parse(Float64,p[2]))
    end
    return qs, ns
end

# x 轴直接使用 k 标签 0..14
xs(ks) = ks .+ 1    # 1-indexed for Plots

# ============================================================
# ── Np=12 能谱对比 ──
# ============================================================
ks_p12, es_p12 = read_spectrum("spectrum_proj_tilted30_Np12.dat")
ks_f12, es_f12 = read_spectrum("../../case3_30sites_Np12/output_Np12/spectrum_Np12.dat")

E0f12 = 0.0   # full ED 已存 E-E0
xs_p12 = xs(ks_p12);  xs_f12 = xs(ks_f12)

# 颜色：最低 Nk=15 个态标红（近15重简并基态流形）
# Np=12 全ED：15态跨 ΔE=0~0.034，之后 gap 到 0.59
# Np=12 投影ED：15态跨 ΔE=0~0.51，之后 gap 到 0.54
gs_idx_f12 = Set(sortperm(es_f12)[1:Nk])
gs_idx_p12 = Set(sortperm(es_p12)[1:Nk])
c_p12 = [i in gs_idx_p12 ? :red       : :steelblue for i in 1:length(es_p12)]
c_f12 = [i in gs_idx_f12 ? :firebrick : :gray50    for i in 1:length(es_f12)]
sz_p12 = [i in gs_idx_p12 ? 7 : 4 for i in 1:length(es_p12)]
sz_f12 = [i in gs_idx_f12 ? 7 : 4 for i in 1:length(es_f12)]

# y 轴范围：显示基态流形 + 第一激发带底部
gap12_f = sort(es_f12)[Nk+1]   # 第16个能级（激发带底）
gap12_p = sort(es_p12)[Nk+1]
emax12  = max(gap12_f, gap12_p) * 1.25

p12 = scatter(xs_f12, es_f12;
    marker=:circle, markersize=sz_f12, markerstrokewidth=0,
    color=c_f12, label="Full 2-band ED (V1=10)",
    xlabel="momentum k",
    ylabel="E − E₀",
    title="TiltedLat30, Np=12 (ν=4/5), 近15重简并基态流形（红色）\nt=1, t'=0.2, V1=10, V2=2, V3=2",
    xticks=(1:Nk, klabels), xrotation=0,
    ylims=(-0.02, emax12), xlims=(0.3, Nk+0.7),
    legend=:topright, framestyle=:box,
    size=(640, 520),
    left_margin=8Plots.mm, bottom_margin=10Plots.mm
)
scatter!(p12, xs_p12, es_p12;
    marker=:diamond, markersize=sz_p12, markerstrokewidth=0,
    color=c_p12, label="Projected ED (lowest band)"
)
# 标注流形顶部与gap
e_top_f12 = sort(es_f12)[Nk]
e_top_p12 = sort(es_p12)[Nk]
annotate!(p12, Nk*0.55, e_top_f12 + emax12*0.06,
    text(@sprintf("Full GS spread=%.4f\ngap→%.4f", e_top_f12, gap12_f), 8, :firebrick, :left))
annotate!(p12, Nk*0.55, e_top_p12 - emax12*0.10,
    text(@sprintf("Proj GS spread=%.4f\ngap→%.4f", e_top_p12, gap12_p), 8, :red, :left))
savefig(p12, "spectrum_proj_vs_full_Np12.pdf")
println("保存：spectrum_proj_vs_full_Np12.pdf")

# ============================================================
# ── Np=13 能谱对比 ──
# ============================================================
ks_p13, es_p13 = read_spectrum("spectrum_proj_tilted30_Np13.dat")
ks_f13, es_f13 = read_spectrum("../../case3_30sites_Np12/output_Np13/spectrum_Np13.dat")

xs_p13 = xs(ks_p13);  xs_f13 = xs(ks_f13)

gs_idx_f13 = Set(sortperm(es_f13)[1:Nk])
gs_idx_p13 = Set(sortperm(es_p13)[1:Nk])
c_p13 = [i in gs_idx_p13 ? :red       : :steelblue for i in 1:length(es_p13)]
c_f13 = [i in gs_idx_f13 ? :firebrick : :gray50    for i in 1:length(es_f13)]
sz_p13 = [i in gs_idx_p13 ? 7 : 4 for i in 1:length(es_p13)]
sz_f13 = [i in gs_idx_f13 ? 7 : 4 for i in 1:length(es_f13)]

gap13_f = sort(es_f13)[Nk+1]
gap13_p = sort(es_p13)[Nk+1]
emax13  = max(gap13_f, gap13_p) * 1.25

p13 = scatter(xs_f13, es_f13;
    marker=:circle, markersize=sz_f13, markerstrokewidth=0,
    color=c_f13, label="Full 2-band ED (V1=10)",
    xlabel="momentum k",
    ylabel="E − E₀",
    title="TiltedLat30, Np=13 (ν=13/15), 近15重简并基态流形（红色）\nt=1, t'=0.2, V1=10, V2=2, V3=2",
    xticks=(1:Nk, klabels), xrotation=0,
    ylims=(-0.02, emax13), xlims=(0.3, Nk+0.7),
    legend=:topright, framestyle=:box,
    size=(640, 520),
    left_margin=8Plots.mm, bottom_margin=10Plots.mm
)
scatter!(p13, xs_p13, es_p13;
    marker=:diamond, markersize=sz_p13, markerstrokewidth=0,
    color=c_p13, label="Projected ED (lowest band)"
)
e_top_f13 = sort(es_f13)[Nk]
e_top_p13 = sort(es_p13)[Nk]
annotate!(p13, Nk*0.55, e_top_f13 + emax13*0.06,
    text(@sprintf("Full GS spread=%.4f\ngap→%.4f", e_top_f13, gap13_f), 8, :firebrick, :left))
annotate!(p13, Nk*0.55, e_top_p13 - emax13*0.10,
    text(@sprintf("Proj GS spread=%.4f\ngap→%.4f", e_top_p13, gap13_p), 8, :red, :left))
savefig(p13, "spectrum_proj_vs_full_Np13.pdf")
println("保存：spectrum_proj_vs_full_Np13.pdf")

# ============================================================
# ── 结构因子 N(q) 对比 ──
# ============================================================
function sq_compare_plot(Np, sq_proj_file, sq_full_file, title_str, outfile)
    qs_p, ns_p = read_sq(sq_proj_file)
    qs_f, ns_f = read_sq(sq_full_file)

    ymax = max(maximum(ns_p), maximum(ns_f)) * 1.20

    p = plot(qs_f .+ 1, ns_f;
        marker=:circle, markersize=7, markerstrokewidth=0.5,
        color=:gray40, linewidth=1.5, label="Full 2-band ED",
        xlabel="momentum q",
        ylabel="N(q)",
        title=title_str,
        xticks=(1:Nk, klabels),
        ylims=(-0.02, ymax), xlims=(0.3, Nk+0.7),
        legend=:topright, framestyle=:box,
        size=(640, 520),
        left_margin=8Plots.mm, bottom_margin=10Plots.mm
    )
    plot!(p, qs_p .+ 1, ns_p;
        marker=:diamond, markersize=7, markerstrokewidth=0.5,
        color=:steelblue, linewidth=1.5, label="Projected ED (lowest band)"
    )
    # 标注两者各自最大峰
    top_f = argmax(ns_f)
    top_p = argmax(ns_p)
    annotate!(p, qs_f[top_f]+1, ns_f[top_f]+ymax*0.05,
        text(@sprintf("%.3f\n(q=%d)", ns_f[top_f], qs_f[top_f]), 8, :gray20, :center))
    annotate!(p, qs_p[top_p]+1+0.5, ns_p[top_p]+ymax*0.05,
        text(@sprintf("%.3f\n(q=%d)", ns_p[top_p], qs_p[top_p]), 8, :steelblue, :center))
    savefig(p, outfile)
    println("保存：$outfile")
end

sq_compare_plot(12,
    "sq_proj_tilted30_Np12.dat",
    "../../case3_30sites_Np12/output_Np12/sq_Np12.dat",
    "Structure Factor N(q) — TiltedLat30, Np=12 (ν=4/5)\nt=1, t'=0.2, V1=10, V2=2, V3=2\n(Full ED kpoints 为旧版本，标签对应，物理 q 略有差异)",
    "sq_proj_vs_full_Np12.pdf"
)

sq_compare_plot(13,
    "sq_proj_tilted30_Np13.dat",
    "../../case3_30sites_Np12/output_Np13/sq_Np13.dat",
    "Structure Factor N(q) — TiltedLat30, Np=13 (ν=13/15)\nt=1, t'=0.2, V1=10, V2=2, V3=2\n(Full ED kpoints 为旧版本，标签对应，物理 q 略有差异)",
    "sq_proj_vs_full_Np13.pdf"
)

# ============================================================
# ── 全扇区结构因子热图 ──
# ============================================================
function sq_all_heatmap(dat_file, title_str, outfile)
    sks=Int[]; qks=Int[]; vals=Float64[]
    for l in readlines(dat_file)
        startswith(l,"#") && continue
        p=split(strip(l)); length(p)<3 && continue
        push!(sks,parse(Int,p[1])); push!(qks,parse(Int,p[2]))
        push!(vals,parse(Float64,p[3]))
    end
    Z = fill(NaN, Nk, Nk)
    for (sk,qk,v) in zip(sks,qks,vals)
        Z[sk+1, qk+1] = v
    end
    ph = heatmap(1:Nk, 1:Nk, Z';
        xlabel="sector k", ylabel="q",
        title=title_str,
        xticks=(1:Nk, klabels), yticks=(1:Nk, klabels),
        color=:viridis, clims=(0, maximum(filter(!isnan,Z))),
        aspect_ratio=:equal, framestyle=:box,
        size=(600, 560),
        left_margin=10Plots.mm, bottom_margin=10Plots.mm
    )
    savefig(ph, outfile)
    println("保存：$outfile")
end

sq_all_heatmap("sq_all_proj_tilted30_Np12.dat",
    "Projected ED N(q) all sectors — TiltedLat30 Np=12",
    "sq_all_proj_Np12.pdf")

sq_all_heatmap("sq_all_proj_tilted30_Np13.dat",
    "Projected ED N(q) all sectors — TiltedLat30 Np=13",
    "sq_all_proj_Np13.pdf")

println("\n完成！共生成 6 个 PDF：")
println("  spectrum_proj_vs_full_Np12.pdf   spectrum_proj_vs_full_Np13.pdf")
println("  sq_proj_vs_full_Np12.pdf         sq_proj_vs_full_Np13.pdf")
println("  sq_all_proj_Np12.pdf             sq_all_proj_Np13.pdf")
