# ============================================================
# 投影 ED 能谱图 — RectLat4x6, Np=4
# 对比 case4 (V1=1) 和 case5 (V1=0) 的 projected ED vs 全 ED
# 用法：julia plot.jl（需先跑完 validate_4x6.jl）
# ============================================================
using Plots, Printf
gr()

Np  = 4
Nk  = 12   # RectLat4x6
n_gs = Nk ÷ Np   # 1/3 filling → 3-fold GS degeneracy

# ── 工具函数 ──
function read_spectrum(path)
    ks = Int[]; es = Float64[]
    for line in readlines(path)
        startswith(line, "#") && continue
        parts = split(strip(line))
        length(parts) < 2 && continue
        push!(ks, parse(Int, parts[1]))
        push!(es, parse(Float64, parts[2]))
    end
    return ks, es
end

# x 轴：(m1, m2) with m1=k%4, m2=k÷4，按 m1 外层 m2 内层排列
order   = [(m1, m2) for m1 in 0:3 for m2 in 0:2]
klabels = ["($m1,$m2)" for (m1, m2) in order]
new_x   = Dict((m1, m2) => i for (i, (m1, m2)) in enumerate(order))
kx(k)   = new_x[(k % 4, k ÷ 4)]

# ── Case 4: V1=1 ──
ks_p4, es_p4 = read_spectrum("spectrum_proj_4x6_V1.dat")
ks_f4, es_f4 = read_spectrum("../../cases/case4_4x6_Np4/spectrum_Np4.dat")

xs_p4 = kx.(ks_p4);  xs_f4 = kx.(ks_f4)

gs_p4 = Set(sortperm(es_p4)[1:n_gs])
gs_f4 = Set(sortperm(es_f4)[1:n_gs])
c_p4  = [i in gs_p4 ? :red    : :steelblue for i in 1:length(es_p4)]
c_f4  = [i in gs_f4 ? :firebrick : :gray60 for i in 1:length(es_f4)]
sz_p4 = [i in gs_p4 ? 9 : 5 for i in 1:length(es_p4)]
sz_f4 = [i in gs_f4 ? 9 : 5 for i in 1:length(es_f4)]

emax4 = 0.35
p_comp4 = scatter(xs_f4, es_f4;
    marker=:circle, markersize=sz_f4, markerstrokewidth=0,
    color=c_f4, label="Full 2-band ED",
    xlabel="momentum (m₁, m₂)",
    ylabel="E − E₀",
    title="RectLat4×6, Np=$Np (ν=1/3),  V₁=1, V₂=V₃=0,  t'=0.2\nProjected ED (blue/red) vs Full 2-band ED (gray/dark)",
    xticks=(1:Nk, klabels), xrotation=45,
    ylims=(-0.01, emax4), xlims=(0.3, Nk+0.7),
    legend=:topright, framestyle=:box,
    size=(640, 520),
    left_margin=8Plots.mm, bottom_margin=12Plots.mm
)
scatter!(p_comp4, xs_p4, es_p4;
    marker=:diamond, markersize=sz_p4, markerstrokewidth=0,
    color=c_p4, label="Projected ED (lowest band)"
)
# 标注投影 ED 前3态
for i in sortperm(es_p4)[1:n_gs]
    annotate!(p_comp4, xs_p4[i], es_p4[i] + 0.010,
        text(@sprintf("%.5f", es_p4[i]), 6, :red, :center))
end
# 标注全 ED 前3态
for i in sortperm(es_f4)[1:n_gs]
    annotate!(p_comp4, xs_f4[i] + 0.35, es_f4[i] + 0.005,
        text(@sprintf("%.5f", es_f4[i]), 6, :firebrick, :center))
end

savefig(p_comp4, "spectrum_proj_vs_full_V1.pdf")
println("保存：spectrum_proj_vs_full_V1.pdf")

# ── Case 4: 纯投影 ED 能谱（单独版）──
p_proj4 = scatter(xs_p4, es_p4;
    marker=:circle, markersize=sz_p4, markerstrokewidth=0,
    color=c_p4,
    xlabel="momentum (m₁, m₂)",
    ylabel="E − E₀",
    title="Projected ED — RectLat4×6, Np=$Np (ν=1/3)\nt=1, t'=0.2, V₁=1, V₂=V₃=0",
    xticks=(1:Nk, klabels), xrotation=45,
    ylims=(-0.01, emax4), xlims=(0.3, Nk+0.7),
    legend=false, framestyle=:box,
    size=(540, 520),
    left_margin=8Plots.mm, bottom_margin=12Plots.mm
)
for i in sortperm(es_p4)[1:n_gs]
    annotate!(p_proj4, xs_p4[i], es_p4[i] + 0.011,
        text(@sprintf("%.5f", es_p4[i]), 7, :red, :center))
end
# 打印 gap 信息
es_sorted = sort(es_p4)
spread = es_sorted[n_gs] - es_sorted[1]
gap34  = es_sorted[n_gs+1] - es_sorted[n_gs]
annotate!(p_proj4, 6.5, emax4*0.88,
    text(@sprintf("spread=%.5f\ngap=%.5f", spread, gap34), 8, :black, :center))
savefig(p_proj4, "spectrum_proj_V1.pdf")
println("保存：spectrum_proj_V1.pdf")

# ── Case 5: V1=0 对比（验证单粒子极限）──
ks_p5, es_p5 = read_spectrum("spectrum_proj_4x6_V0.dat")
ks_f5, es_f5 = read_spectrum("../../cases/case5_4x6_Np4_V0/spectrum_Np4.dat")
xs_p5 = kx.(ks_p5);  xs_f5 = kx.(ks_f5)

N_show = min(length(es_p5), length(es_f5), 48)
emax5  = max(es_p5[N_show], es_f5[N_show]) * 1.15

p_comp5 = scatter(xs_f5[1:N_show], es_f5[1:N_show];
    marker=:circle, markersize=5, markerstrokewidth=0,
    color=:gray50, label="Full 2-band ED",
    xlabel="momentum (m₁, m₂)",
    ylabel="E − E₀",
    title="RectLat4×6, Np=$Np,  V₁=0 (non-interacting)\nProjected ED vs Full 2-band ED",
    xticks=(1:Nk, klabels), xrotation=45,
    ylims=(-0.01, emax5), xlims=(0.3, Nk+0.7),
    legend=:topright, framestyle=:box,
    size=(640, 520),
    left_margin=8Plots.mm, bottom_margin=12Plots.mm
)
scatter!(p_comp5, xs_p5[1:N_show], es_p5[1:N_show];
    marker=:diamond, markersize=4, markerstrokewidth=0,
    color=:steelblue, label="Projected ED (lowest band)"
)
savefig(p_comp5, "spectrum_proj_vs_full_V0.pdf")
println("保存：spectrum_proj_vs_full_V0.pdf")

# ============================================================
# 结构因子 N(q) 对比图
# ============================================================

function read_sq(path)
    qs = Int[]; ns = Float64[]
    for line in readlines(path)
        startswith(line, "#") && continue
        parts = split(strip(line))
        length(parts) < 2 && continue
        push!(qs, parse(Int, parts[1]))
        push!(ns, parse(Float64, parts[2]))
    end
    return qs, ns
end

function read_sq_all(path)
    sks = Int[]; qks = Int[]; vals = Float64[]
    for line in readlines(path)
        startswith(line, "#") && continue
        parts = split(strip(line))
        length(parts) < 3 && continue
        push!(sks, parse(Int, parts[1]))
        push!(qks, parse(Int, parts[2]))
        push!(vals, parse(Float64, parts[3]))
    end
    return sks, qks, vals
end

# ── Case 4: N(q) 投影 ED vs 全 ED ──
qs_p4, ns_p4 = read_sq("sq_proj_4x6_V1.dat")
qs_f4, ns_f4 = read_sq("../../cases/case4_4x6_Np4/sq_Np4.dat")

xs_sq_p4 = kx.(qs_p4)
xs_sq_f4 = kx.(qs_f4)
ord_p4   = sortperm(xs_sq_p4)
ord_f4   = sortperm(xs_sq_f4)

sq_ymax = max(maximum(ns_f4), maximum(ns_p4)) * 1.20

p_sq4 = plot(xs_sq_f4[ord_f4], ns_f4[ord_f4];
    marker=:circle, markersize=7, markerstrokewidth=0.5,
    color=:gray40, linewidth=1.5, label="Full 2-band ED",
    xlabel="momentum q  (m₁, m₂)",
    ylabel="N(q)",
    title="Structure Factor N(q) — RectLat4×6, Np=$Np (ν=1/3)\nt=1, t'=0.2, V₁=1",
    xticks=(1:Nk, klabels), xrotation=45,
    ylims=(-0.02, sq_ymax), xlims=(0.3, Nk+0.7),
    legend=:topright, framestyle=:box,
    size=(640, 520),
    left_margin=8Plots.mm, bottom_margin=14Plots.mm
)
plot!(p_sq4, xs_sq_p4[ord_p4], ns_p4[ord_p4];
    marker=:diamond, markersize=7, markerstrokewidth=0.5,
    color=:steelblue, linewidth=1.5, label="Projected ED (lowest band)"
)
# 标注全 ED 最大两个峰
top2_f4 = sort(findall(q -> q != 0, qs_f4), by=i -> -ns_f4[i])[1:2]
for i in top2_f4
    annotate!(p_sq4, xs_sq_f4[i], ns_f4[i] + sq_ymax*0.04,
        text(@sprintf("%.4f\n(k=%d)", ns_f4[i], qs_f4[i]), 7, :gray20, :center))
end
# 标注投影 ED 最大两个峰
top2_p4 = sort(findall(q -> q != 0, qs_p4), by=i -> -ns_p4[i])[1:2]
for i in top2_p4
    annotate!(p_sq4, xs_sq_p4[i] + 0.45, ns_p4[i] + sq_ymax*0.04,
        text(@sprintf("%.4f\n(q=%d)", ns_p4[i], qs_p4[i]), 7, :steelblue, :center))
end
savefig(p_sq4, "sq_proj_vs_full_V1.pdf")
println("保存：sq_proj_vs_full_V1.pdf")

# ── Case 5: N(q) V1=0 验证 ──
qs_p5, ns_p5 = read_sq("sq_proj_4x6_V0.dat")
qs_f5, ns_f5 = read_sq("../../cases/case5_4x6_Np4_V0/sq_Np4.dat")

xs_sq_p5 = kx.(qs_p5)
xs_sq_f5 = kx.(qs_f5)
ord_p5   = sortperm(xs_sq_p5)
ord_f5   = sortperm(xs_sq_f5)

sq_ymax5 = max(maximum(ns_f5), maximum(ns_p5)) * 1.20

p_sq5 = plot(xs_sq_f5[ord_f5], ns_f5[ord_f5];
    marker=:circle, markersize=7, markerstrokewidth=0.5,
    color=:gray40, linewidth=1.5, label="Full 2-band ED",
    xlabel="momentum q  (m₁, m₂)",
    ylabel="N(q)",
    title="Structure Factor N(q) — RectLat4×6, Np=$Np,  V₁=0 (non-interacting)\nt=1, t'=0.2",
    xticks=(1:Nk, klabels), xrotation=45,
    ylims=(-0.005, sq_ymax5), xlims=(0.3, Nk+0.7),
    legend=:topright, framestyle=:box,
    size=(640, 520),
    left_margin=8Plots.mm, bottom_margin=14Plots.mm
)
plot!(p_sq5, xs_sq_p5[ord_p5], ns_p5[ord_p5];
    marker=:diamond, markersize=7, markerstrokewidth=0.5,
    color=:steelblue, linewidth=1.5, label="Projected ED (lowest band)"
)
savefig(p_sq5, "sq_proj_vs_full_V0.pdf")
println("保存：sq_proj_vs_full_V0.pdf")

# ── 全扇区热图：投影 ED  vs 全 ED ──
sks_p, qks_p, vals_p = read_sq_all("sq_all_proj_4x6_V1.dat")
sks_f, qks_f, vals_f = read_sq_all("../../cases/case4_4x6_Np4/sq_all_Np4.dat")

Zp = fill(NaN, Nk, Nk)
for (sk, qk, v) in zip(sks_p, qks_p, vals_p)
    Zp[kx(sk), kx(qk)] = v
end

Zf = fill(NaN, Nk, Nk)
for (sk, qk, v) in zip(sks_f, qks_f, vals_f)
    Zf[kx(sk), kx(qk)] = v
end

clim_max = max(maximum(filter(!isnan, Zp)), maximum(filter(!isnan, Zf)))

ph1 = heatmap(1:Nk, 1:Nk, Zp';
    xlabel="sector k", ylabel="q",
    title="Projected ED — N(q) all sectors\nV₁=1",
    xticks=(1:Nk, klabels), yticks=(1:Nk, klabels),
    xrotation=45, color=:viridis,
    clims=(0, clim_max),
    aspect_ratio=:equal, framestyle=:box,
    size=(600, 560),
    left_margin=14Plots.mm, bottom_margin=14Plots.mm
)

ph2 = heatmap(1:Nk, 1:Nk, Zf';
    xlabel="sector k", ylabel="q",
    title="Full 2-band ED — N(q) all sectors\nV₁=1",
    xticks=(1:Nk, klabels), yticks=(1:Nk, klabels),
    xrotation=45, color=:viridis,
    clims=(0, clim_max),
    aspect_ratio=:equal, framestyle=:box,
    size=(600, 560),
    left_margin=14Plots.mm, bottom_margin=14Plots.mm
)

savefig(ph1, "sq_all_proj_V1.pdf")
savefig(ph2, "sq_all_full_V1.pdf")
println("保存：sq_all_proj_V1.pdf")
println("保存：sq_all_full_V1.pdf")

println("\n完成！输出文件：")
println("  spectrum_proj_vs_full_V1.pdf  — Case4 投影 vs 全 ED 能谱对比")
println("  spectrum_proj_V1.pdf          — Case4 投影 ED 能谱单独版")
println("  spectrum_proj_vs_full_V0.pdf  — Case5 V1=0 验证")
println("  sq_proj_vs_full_V1.pdf        — Case4 N(q) 投影 vs 全 ED 对比")
println("  sq_proj_vs_full_V0.pdf        — Case5 N(q) V1=0 对比")
println("  sq_all_proj_V1.pdf            — Case4 全扇区 N(q) 热图（投影）")
println("  sq_all_full_V1.pdf            — Case4 全扇区 N(q) 热图（全 ED）")
