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

# 颜色：基态标红
n_gs_p12 = sum(e < 1e-8 for e in es_p12)
n_gs_f12 = sum(e < 1e-8 for e in es_f12)
c_p12 = [e < 1e-8 ? :red       : :steelblue for e in es_p12]
c_f12 = [e < 1e-8 ? :firebrick : :gray50    for e in es_f12]

emax12 = 0.8

p12 = scatter(xs_f12, es_f12;
    marker=:circle, markersize=6, markerstrokewidth=0,
    color=c_f12, label="Full 2-band ED (V1=10)",
    xlabel="momentum k",
    ylabel="E − E₀",
    title="TiltedLat30, Np=12 (ν=4/5)\nt=1, t'=0.2, V1=10, V2=2, V3=2",
    xticks=(1:Nk, klabels), xrotation=0,
    ylims=(-0.02, emax12), xlims=(0.3, Nk+0.7),
    legend=:topright, framestyle=:box,
    size=(640, 520),
    left_margin=8Plots.mm, bottom_margin=10Plots.mm
)
scatter!(p12, xs_p12, es_p12;
    marker=:diamond, markersize=6, markerstrokewidth=0,
    color=c_p12, label="Projected ED (lowest band)"
)
# 标注基态扇区
for (m,e) in zip(ks_f12, es_f12)
    e < 1e-8 || continue
    annotate!(p12, m+1, e+0.025, text("k=$m", 8, :firebrick, :center))
end
for (m,e) in zip(ks_p12, es_p12)
    e < 1e-8 || continue
    annotate!(p12, m+1+0.4, e+0.012, text("k=$m", 7, :red, :center))
end
savefig(p12, "spectrum_proj_vs_full_Np12.pdf")
println("保存：spectrum_proj_vs_full_Np12.pdf")

# ============================================================
# ── Np=13 能谱对比 ──
# ============================================================
ks_p13, es_p13 = read_spectrum("spectrum_proj_tilted30_Np13.dat")
ks_f13, es_f13 = read_spectrum("../../case3_30sites_Np12/output_Np13/spectrum_Np13.dat")

xs_p13 = xs(ks_p13);  xs_f13 = xs(ks_f13)

c_p13 = [e < 1e-8 ? :red       : :steelblue for e in es_p13]
c_f13 = [e < 1e-8 ? :firebrick : :gray50    for e in es_f13]

emax13 = min(0.6, max(maximum(es_p13[1:min(30,end)]), maximum(es_f13[1:min(30,end)]))*1.2)

p13 = scatter(xs_f13, es_f13;
    marker=:circle, markersize=6, markerstrokewidth=0,
    color=c_f13, label="Full 2-band ED (V1=10)",
    xlabel="momentum k",
    ylabel="E − E₀",
    title="TiltedLat30, Np=13 (ν=13/15)\nt=1, t'=0.2, V1=10, V2=2, V3=2",
    xticks=(1:Nk, klabels), xrotation=0,
    ylims=(-0.02, emax13), xlims=(0.3, Nk+0.7),
    legend=:topright, framestyle=:box,
    size=(640, 520),
    left_margin=8Plots.mm, bottom_margin=10Plots.mm
)
scatter!(p13, xs_p13, es_p13;
    marker=:diamond, markersize=6, markerstrokewidth=0,
    color=c_p13, label="Projected ED (lowest band)"
)
for (m,e) in zip(ks_f13, es_f13)
    e < 1e-8 || continue
    annotate!(p13, m+1, e+emax13*0.07, text("k=$m", 8, :firebrick, :center))
end
for (m,e) in zip(ks_p13, es_p13)
    e < 1e-8 || continue
    annotate!(p13, m+1+0.4, e+emax13*0.04, text("k=$m", 7, :red, :center))
end
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
