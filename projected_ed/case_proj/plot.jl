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
ks_f4, es_f4 = read_spectrum("../../case4_4x6_Np4/spectrum_Np4.dat")

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
ks_f5, es_f5 = read_spectrum("../../case5_4x6_Np4_V0/spectrum_Np4.dat")
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

println("\n完成！输出文件：")
println("  spectrum_proj_vs_full_V1.pdf  — Case4 投影 vs 全 ED 对比")
println("  spectrum_proj_V1.pdf          — Case4 投影 ED 单独版")
println("  spectrum_proj_vs_full_V0.pdf  — Case5 V1=0 验证")
